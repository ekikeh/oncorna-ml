"""Audit TCGA-BRCA identifiers and describe expression/clinical overlap.

This script reads the expression matrix one row at a time. It does not filter,
normalize, merge, or otherwise alter either source file. From the repository
root, run:

    python scripts/audit_data.py

The generated report is written to ``docs/data_audit.md``.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import re
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from zoneinfo import ZoneInfo, ZoneInfoNotFoundError

PROJECT_ROOT = Path(__file__).resolve().parents[1]
EXPRESSION_PATH = PROJECT_ROOT / "data" / "raw" / "TCGA-BRCA_HiSeqV2.tsv.gz"
CLINICAL_PATH = PROJECT_ROOT / "data" / "raw" / "TCGA-BRCA_BRCA_clinicalMatrix.tsv"
REPORT_PATH = PROJECT_ROOT / "docs" / "data_audit.md"
XENA_PAGE_GENE_ID_COUNT = 20_531
XENA_EXPRESSION_METADATA_URL = (
    "https://xenabrowser.net/datapages/?dataset=TCGA.BRCA.sampleMap%2FHiSeqV2"
    "&host=https%3A%2F%2Ftcga.xenahubs.net"
)

MISSING_VALUES = {
    "",
    "--",
    "n/a",
    "na",
    "nan",
    "none",
    "not applicable",
    "not available",
    "null",
}
ENSEMBL_GENE_ID = re.compile(r"^ENSG\d+(?:\.\d+)?$", re.IGNORECASE)
TCGA_CENTER_CODE = re.compile(r"^[A-Z0-9]{2}$", re.IGNORECASE)
TCGA_PARTICIPANT_CODE = re.compile(r"^[A-Z0-9]{4}$", re.IGNORECASE)


@dataclass
class ExpressionAudit:
    identifier_column: str
    sample_ids: list[str]
    sample_columns: int
    empty_sample_columns: int
    duplicate_sample_ids: dict[str, int]
    gene_rows: int
    nonempty_gene_ids: int
    unique_gene_ids: int
    duplicate_gene_ids: dict[str, int]
    duplicate_gene_rows: int
    first_gene_ids: list[str]
    identifier_types: dict[str, int]
    malformed_rows: int
    malformed_row_examples: list[int]


@dataclass
class ClinicalAudit:
    columns: list[str]
    row_count: int
    malformed_rows: int
    malformed_row_examples: list[int]
    sample_id_column: str | None
    sample_ids: list[str]
    missing_sample_ids: int
    duplicate_sample_ids: dict[str, int]
    key_columns: dict[str, list[str]]
    missing_counts: dict[str, int]
    pam50_column: str | None
    pam50_counts: dict[str, int]
    pam50_missing: int
    pam50_case_or_spacing_variants: list[list[str]]
    pam50_other_field_counts: dict[str, dict[str, int]]
    patient_suffix_comparisons: int
    patient_suffix_matches: int
    patient_suffix_examples: list[tuple[str, str, str]]


@dataclass
class SampleMatchAudit:
    expression_unique_ids: int
    clinical_unique_ids: int
    intersection_ids: list[str]
    expression_only_ids: list[str]
    clinical_only_ids: list[str]


@dataclass
class PatientAudit:
    patient_sample_counts: dict[str, int]
    patients_with_multiple_samples: int
    samples_per_patient_distribution: dict[int, int]
    duplicate_patient_examples: list[tuple[str, list[str]]]
    unparseable_sample_ids: list[str]


@dataclass
class AuditResult:
    expression: ExpressionAudit
    clinical: ClinicalAudit
    matching: SampleMatchAudit
    patients: PatientAudit
    markdown: str


def is_missing(value: str | None) -> bool:
    """Return True for explicit blank/null-style values, not biological categories."""
    return value is None or value.strip().casefold() in MISSING_VALUES


def parse_tcga_patient_id(sample_barcode: str) -> str | None:
    """Take the first three hyphen-separated fields of a TCGA barcode.

    For example, ``TCGA-AR-A5QQ-01`` becomes ``TCGA-AR-A5QQ``. The returned
    identifier is the participant-level prefix; sample/aliquot suffixes are
    not included. Original barcode case is preserved.
    """
    parts = sample_barcode.strip().split("-")
    if len(parts) < 3 or parts[0].upper() != "TCGA":
        return None
    if not TCGA_CENTER_CODE.fullmatch(parts[1]):
        return None
    if not TCGA_PARTICIPANT_CODE.fullmatch(parts[2]):
        return None
    return "-".join(parts[:3])


def classify_gene_identifier(identifier: str) -> str:
    """Describe the visible gene-ID pattern without remapping any identifiers."""
    value = identifier.strip()
    if not value or value.casefold() in MISSING_VALUES:
        return "Blank or missing identifier"

    candidate = value
    is_xena_fallback = value.startswith("?|")
    if is_xena_fallback:
        candidate = value.split("|", maxsplit=1)[1]

    if ENSEMBL_GENE_ID.fullmatch(candidate):
        if is_xena_fallback:
            return "Xena fallback containing an Ensembl ID"
        return "Ensembl gene ID"
    if is_xena_fallback and candidate.isdigit():
        return "Xena fallback with numeric ID"
    if value.startswith("?"):
        return "Other Xena/unmapped identifier"
    if candidate.isdigit():
        return "Numeric identifier"
    return "Gene-symbol-like identifier"


def _duplicate_values(values: list[str]) -> dict[str, int]:
    counts = Counter(value for value in values if value)
    return dict(sorted((value, count) for value, count in counts.items() if count > 1))


def inspect_expression(path: Path) -> ExpressionAudit:
    """Stream the gzipped expression matrix and summarize IDs and dimensions."""
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(
            f"Expression input not found: {path}. Run scripts/download_data.py first."
        )

    try:
        with gzip.open(path, "rt", encoding="utf-8-sig", newline="") as handle:
            reader = csv.reader(handle, delimiter="\t")
            header = next(reader, None)
            if not header or len(header) < 2:
                raise ValueError(f"Expression file has no usable header: {path}")

            sample_ids = [sample.strip() for sample in header[1:] if sample.strip()]
            empty_sample_columns = sum(not sample.strip() for sample in header[1:])
            gene_counts: Counter[str] = Counter()
            identifier_types: Counter[str] = Counter()
            first_gene_ids: list[str] = []
            gene_rows = 0
            nonempty_gene_ids = 0
            malformed_rows = 0
            malformed_row_examples: list[int] = []

            for row_number, row in enumerate(reader, start=2):
                gene_rows += 1
                if len(row) != len(header):
                    malformed_rows += 1
                    if len(malformed_row_examples) < 10:
                        malformed_row_examples.append(row_number)
                if not row:
                    continue
                gene_id = row[0].strip()
                if not gene_id:
                    continue
                nonempty_gene_ids += 1
                gene_counts[gene_id] += 1
                identifier_types[classify_gene_identifier(gene_id)] += 1
                if len(first_gene_ids) < 10:
                    first_gene_ids.append(gene_id)
    except (OSError, EOFError, csv.Error) as error:
        raise ValueError(f"Could not read expression matrix {path}: {error}") from error

    duplicate_gene_ids = {
        gene_id: count for gene_id, count in sorted(gene_counts.items()) if count > 1
    }
    duplicate_sample_ids = _duplicate_values(sample_ids)
    return ExpressionAudit(
        identifier_column=header[0].strip(),
        sample_ids=sample_ids,
        sample_columns=len(header) - 1,
        empty_sample_columns=empty_sample_columns,
        duplicate_sample_ids=duplicate_sample_ids,
        gene_rows=gene_rows,
        nonempty_gene_ids=nonempty_gene_ids,
        unique_gene_ids=len(gene_counts),
        duplicate_gene_ids=duplicate_gene_ids,
        duplicate_gene_rows=sum(count - 1 for count in duplicate_gene_ids.values()),
        first_gene_ids=first_gene_ids,
        identifier_types=dict(sorted(identifier_types.items())),
        malformed_rows=malformed_rows,
        malformed_row_examples=malformed_row_examples,
    )


def _key_columns(columns: list[str]) -> dict[str, list[str]]:
    """Find likely key fields by name; preserve the source's exact spelling."""
    categories: dict[str, list[str]] = {
        "Sample ID": [],
        "Patient ID": [],
        "PAM50 subtype": [],
        "ER": [],
        "PR": [],
        "HER2": [],
        "Age": [],
        "Race/ethnicity": [],
        "Survival": [],
    }

    for column in columns:
        name = column.casefold()
        if name in {"sampleid", "sample_id", "sample", "bcr_sample_barcode"}:
            categories["Sample ID"].append(column)
        if name in {"patient_id", "bcr_patient_barcode", "_patient"}:
            categories["Patient ID"].append(column)
        if "pam50" in name:
            categories["PAM50 subtype"].append(column)
        has_er_status_name = re.search(
            r"(?:^|_)(?:er_status|estrogen_receptor_status)(?:_|$)", name
        )
        has_er_measurement_name = any(
            marker in name
            for marker in (
                "estrogen_receptor",
                "er_detection",
                "er_level",
                "_er_pos_",
                "_er_pos_cell_",
            )
        )
        if has_er_status_name or has_er_measurement_name:
            categories["ER"].append(column)
        if any(
            marker in name
            for marker in (
                "pr_status",
                "progesterone_receptor",
                "prgstrn",
                "pgr_detection",
                "_pr_pos_",
            )
        ):
            categories["PR"].append(column)
        if "her2" in name or "erbb2" in name:
            categories["HER2"].append(column)
        has_age_token = re.search(r"(?:^|_)age(?:_|$)", name)
        if has_age_token or name in {"days_to_birth", "year_of_initial_pathologic_diagnosis"}:
            categories["Age"].append(column)
        if "race" in name or "ethnic" in name:
            categories["Race/ethnicity"].append(column)
        if re.search(r"(?:^|_)(?:os|rfs|dss|dfi|pfs)(?:_|$)", name) or any(
            marker in name
            for marker in (
                "survival",
                "vital_status",
                "days_to_death",
                "days_to_date_of_death",
                "days_to_last",
                "last_contact",
            )
        ):
            categories["Survival"].append(column)

    return categories


def _find_column(columns: list[str], preferred_names: tuple[str, ...]) -> str | None:
    by_casefold = {column.casefold(): column for column in columns}
    for preferred in preferred_names:
        if preferred.casefold() in by_casefold:
            return by_casefold[preferred.casefold()]
    return None


def _label_variant_groups(labels: list[str]) -> list[list[str]]:
    groups: dict[str, set[str]] = defaultdict(set)
    for label in labels:
        groups[" ".join(label.strip().split()).casefold()].add(label)
    return [sorted(group) for group in groups.values() if len(group) > 1]


def inspect_clinical(path: Path) -> ClinicalAudit:
    """Stream the clinical matrix and summarize IDs, labels, and missingness."""
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(
            f"Clinical input not found: {path}. Run scripts/download_data.py first."
        )

    try:
        with path.open("rt", encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            columns = reader.fieldnames
            if not columns:
                raise ValueError(f"Clinical file has no usable header: {path}")

            key_columns = _key_columns(columns)
            sample_id_column = _find_column(columns, ("sampleID", "sample_id", "sample"))
            if sample_id_column is None:
                sample_id_column = _find_column(columns, ("bcr_sample_barcode",))
            pam50_column = _find_column(columns, ("PAM50Call_RNAseq",))
            pam50_columns = key_columns["PAM50 subtype"]
            tracked_columns = sorted({column for group in key_columns.values() for column in group})

            row_count = 0
            malformed_rows = 0
            malformed_row_examples: list[int] = []
            sample_ids: list[str] = []
            missing_sample_ids = 0
            pam50_missing = 0
            missing_counts: Counter[str] = Counter()
            pam50_distributions: dict[str, Counter[str]] = {
                column: Counter() for column in pam50_columns
            }
            patient_suffix_comparisons = 0
            patient_suffix_matches = 0
            patient_suffix_examples: list[tuple[str, str, str]] = []

            for row_number, row in enumerate(reader, start=2):
                row_count += 1
                if None in row or any(row.get(column) is None for column in columns):
                    malformed_rows += 1
                    if len(malformed_row_examples) < 10:
                        malformed_row_examples.append(row_number)

                for column in tracked_columns:
                    value = row.get(column)
                    if is_missing(value):
                        missing_counts[column] += 1

                sample_id_value = row.get(sample_id_column, "") if sample_id_column else ""
                sample_id = (sample_id_value or "").strip()
                if not sample_id:
                    missing_sample_ids += 1
                else:
                    sample_ids.append(sample_id)

                for column in pam50_columns:
                    value = row.get(column)
                    if is_missing(value):
                        pam50_distributions[column]["<MISSING>"] += 1
                    else:
                        pam50_distributions[column][value] += 1

                if pam50_column:
                    value = row.get(pam50_column)
                    if is_missing(value):
                        pam50_missing += 1

                source_patient_id = (row.get("patient_id") or "").strip()
                full_patient_barcode = (row.get("bcr_patient_barcode") or "").strip()
                if source_patient_id and full_patient_barcode:
                    patient_suffix_comparisons += 1
                    if source_patient_id == full_patient_barcode.split("-")[-1]:
                        patient_suffix_matches += 1
                    if len(patient_suffix_examples) < 3:
                        patient_suffix_examples.append(
                            (sample_id, source_patient_id, full_patient_barcode)
                        )
    except (OSError, csv.Error) as error:
        raise ValueError(f"Could not read clinical matrix {path}: {error}") from error

    duplicate_sample_ids = _duplicate_values(sample_ids)
    pam50_counts: dict[str, int] = {}
    if pam50_column:
        distribution = pam50_distributions[pam50_column]
        pam50_counts = {
            label: count for label, count in sorted(distribution.items()) if label != "<MISSING>"
        }
        pam50_missing = distribution.get("<MISSING>", 0)

    other_pam50_counts = {
        column: dict(sorted(counts.items()))
        for column, counts in pam50_distributions.items()
        if column != pam50_column
    }
    pam50_labels = list(pam50_counts)

    return ClinicalAudit(
        columns=columns,
        row_count=row_count,
        malformed_rows=malformed_rows,
        malformed_row_examples=malformed_row_examples,
        sample_id_column=sample_id_column,
        sample_ids=sample_ids,
        missing_sample_ids=missing_sample_ids,
        duplicate_sample_ids=duplicate_sample_ids,
        key_columns=key_columns,
        missing_counts={column: missing_counts[column] for column in tracked_columns},
        pam50_column=pam50_column,
        pam50_counts=pam50_counts,
        pam50_missing=pam50_missing,
        pam50_case_or_spacing_variants=_label_variant_groups(pam50_labels),
        pam50_other_field_counts=other_pam50_counts,
        patient_suffix_comparisons=patient_suffix_comparisons,
        patient_suffix_matches=patient_suffix_matches,
        patient_suffix_examples=patient_suffix_examples,
    )


def match_samples(expression: ExpressionAudit, clinical: ClinicalAudit) -> SampleMatchAudit:
    """Match exact TCGA sample IDs after trimming outer whitespace only."""
    expression_ids = set(expression.sample_ids)
    clinical_ids = set(clinical.sample_ids)
    return SampleMatchAudit(
        expression_unique_ids=len(expression_ids),
        clinical_unique_ids=len(clinical_ids),
        intersection_ids=sorted(expression_ids & clinical_ids),
        expression_only_ids=sorted(expression_ids - clinical_ids),
        clinical_only_ids=sorted(clinical_ids - expression_ids),
    )


def audit_expression_patients(sample_ids: list[str]) -> PatientAudit:
    """Group expression sample barcodes by their first three TCGA fields."""
    samples_by_patient: dict[str, list[str]] = defaultdict(list)
    unparseable: list[str] = []
    for sample_id in sample_ids:
        patient_id = parse_tcga_patient_id(sample_id)
        if patient_id is None:
            unparseable.append(sample_id)
        else:
            samples_by_patient[patient_id].append(sample_id)

    distribution = Counter(len(samples) for samples in samples_by_patient.values())
    duplicate_examples = [
        (patient_id, samples_by_patient[patient_id])
        for patient_id in sorted(samples_by_patient)
        if len(samples_by_patient[patient_id]) > 1
    ][:10]
    return PatientAudit(
        patient_sample_counts=dict(
            sorted((patient_id, len(samples)) for patient_id, samples in samples_by_patient.items())
        ),
        patients_with_multiple_samples=sum(
            patient_count
            for sample_count, patient_count in distribution.items()
            if sample_count > 1
        ),
        samples_per_patient_distribution=dict(sorted(distribution.items())),
        duplicate_patient_examples=duplicate_examples,
        unparseable_sample_ids=unparseable,
    )


def _format_percentage(numerator: int, denominator: int) -> str:
    if denominator == 0:
        return "n/a"
    return f"{numerator / denominator * 100:.1f}%"


def _code_list(values: list[str]) -> str:
    if not values:
        return "None"
    return ", ".join(f"`{value}`" for value in values)


def _render_missingness(clinical: ClinicalAudit) -> list[str]:
    missing_tokens = [f"`{token or 'blank'}`" for token in sorted(MISSING_VALUES)]
    lines = [
        "Missingness treats blank strings and explicit null-style tokens "
        f"({', '.join(missing_tokens)}) as missing. "
        "Values such as `Indeterminate` and `0` are retained as observed values.",
        "",
        "| Category | Exact source column | Present | Missing |",
        "|---|---|---:|---:|",
    ]
    any_rows = False
    for category, columns in clinical.key_columns.items():
        for column in columns:
            any_rows = True
            missing = clinical.missing_counts.get(column, 0)
            present = clinical.row_count - missing
            lines.append(
                f"| {category} | `{column}` | {present}/{clinical.row_count} "
                f"({_format_percentage(present, clinical.row_count)}) | "
                f"{missing}/{clinical.row_count} "
                f"({_format_percentage(missing, clinical.row_count)}) |"
            )
    if not any_rows:
        lines.append("| — | No matching key fields found | — | — |")
    return lines


def _render_report(
    expression: ExpressionAudit,
    clinical: ClinicalAudit,
    matching: SampleMatchAudit,
    patients: PatientAudit,
    audit_date: str,
) -> str:
    expression_set_size = matching.expression_unique_ids
    clinical_set_size = matching.clinical_unique_ids
    intersection_size = len(matching.intersection_ids)

    lines = [
        "# TCGA-BRCA data audit and identifier harmonization",
        "",
        f"**Generated:** {audit_date}<br>",
        (
            "**Scope:** Phase 1B data audit only. No samples or genes were filtered, "
            "no labels were changed, and no modeling/DE/GSEA was performed."
        ),
        "",
        "## Inputs and dimensions",
        "",
        "| Input | Rows/items | Columns/items | Notes |",
        "|---|---:|---:|---|",
        f"| Expression (`{EXPRESSION_PATH.name}`) | {expression.gene_rows:,} gene rows | "
        f"{expression.sample_columns:,} sample columns | First field: "
        f"`{expression.identifier_column}` |",
        f"| Clinical (`{CLINICAL_PATH.name}`) | {clinical.row_count:,} sample rows | "
        f"{len(clinical.columns):,} fields | Tab-separated; sample ID field: "
        f"`{clinical.sample_id_column or 'not found'}` |",
        "",
        "### Expression matrix structure",
        "",
        f"- Gene rows: **{expression.gene_rows:,}**; "
        f"nonblank identifiers: **{expression.nonempty_gene_ids:,}**; "
        f"unique gene IDs: **{expression.unique_gene_ids:,}**.",
        f"- Expression sample columns: **{expression.sample_columns:,}**; "
        f"unique nonblank sample IDs: **{matching.expression_unique_ids:,}**; "
        f"blank sample headers: **{expression.empty_sample_columns}**.",
        f"- First gene identifiers: {_code_list(expression.first_gene_ids[:5])}.",
        f"- First sample IDs: {_code_list(expression.sample_ids[:5])}.",
        "- Gene identifier pattern counts (classifying visible strings only; no mapping applied):",
    ]
    if expression.identifier_types:
        for identifier_type, count in expression.identifier_types.items():
            lines.append(f"  - {identifier_type}: {count:,}")
    else:
        lines.append("  - No nonblank gene identifiers found.")

    ensembl_types = {"Ensembl gene ID", "Xena fallback containing an Ensembl ID"}
    ensembl_identifier_count = sum(
        count
        for identifier_type, count in expression.identifier_types.items()
        if identifier_type in ensembl_types
    )
    lines.append(
        f"- Ensembl-style identifiers detected by this pattern check: "
        f"**{ensembl_identifier_count:,}**."
    )

    if expression.duplicate_gene_ids:
        examples = list(expression.duplicate_gene_ids.items())[:10]
        formatted = ", ".join(f"`{gene}` × {count}" for gene, count in examples)
        lines.append(
            "- Duplicate gene identifiers: "
            f"**{len(expression.duplicate_gene_ids):,} distinct IDs** "
            f"({expression.duplicate_gene_rows:,} extra rows). Examples: {formatted}."
        )
    else:
        lines.append("- Duplicate gene identifiers: **none detected**.")

    if expression.duplicate_sample_ids:
        examples = list(expression.duplicate_sample_ids.items())[:10]
        formatted = ", ".join(f"`{sample}` × {count}" for sample, count in examples)
        lines.append(
            f"- Duplicate exact expression sample IDs: **{len(expression.duplicate_sample_ids)}**. "
            f"Examples: {formatted}."
        )
    else:
        lines.append("- Duplicate exact expression sample IDs: **none detected**.")
    lines.append(f"- Rows with a width different from the header: **{expression.malformed_rows}**.")
    if expression.malformed_row_examples:
        example_rows = _code_list([str(row) for row in expression.malformed_row_examples])
        lines.append(f"  - First affected file row numbers: {example_rows}.")
    if expression.gene_rows != XENA_PAGE_GENE_ID_COUNT:
        difference = abs(expression.gene_rows - XENA_PAGE_GENE_ID_COUNT)
        lines.append(
            f"- Xena Data Pages metadata lists {XENA_PAGE_GENE_ID_COUNT:,} identifiers, "
            f"while the downloaded file contains {expression.gene_rows:,} data rows "
            f"after its header (difference: {difference:,}). This report uses the scanned "
            f"file row count; the cause of the discrepancy is unresolved. "
            f"[Metadata]({XENA_EXPRESSION_METADATA_URL})."
        )

    lines.extend(
        [
            "",
            "### Identifier systems",
            "",
            "- Expression columns use TCGA sample barcodes (for example, `TCGA-AR-A5QQ-01`).",
            (
                "- Expression row identifiers are predominantly gene-symbol-like; any Ensembl-like "
                "or Xena fallback IDs are counted above. This is a pattern audit, not annotation "
                "remapping."
            ),
            (
                "- The clinical matrix exposes `sampleID`, `patient_id`, `bcr_patient_barcode`, "
                "and `_PATIENT`. The clinical `patient_id` field is not the full TCGA participant "
                "barcode: where both fields are present, it matches the final barcode segment. "
                "The full barcode fields retain the `TCGA-XX-XXXX` form."
            ),
            (
                "- For expression barcodes, the participant-level ID uses the first three "
                "hyphen-separated fields joined unchanged: `TCGA-AR-A5QQ-01` → `TCGA-AR-A5QQ`. "
                "This keeps sample/aliquot suffixes out of the patient ID. No samples are "
                "collapsed by this audit."
            ),
            (
                "- Comparison check for clinical `patient_id` suffix versus `bcr_patient_barcode`: "
                f"{clinical.patient_suffix_matches:,} matches among "
                f"{clinical.patient_suffix_comparisons:,} rows with both fields present."
            ),
            (
                "- Example clinical mapping: "
                f"`sampleID` `{clinical.patient_suffix_examples[0][0]}`, "
                f"`patient_id` `{clinical.patient_suffix_examples[0][1]}` matches the final "
                f"segment of `bcr_patient_barcode` `{clinical.patient_suffix_examples[0][2]}`."
                if clinical.patient_suffix_examples
                else "- No rows had both clinical patient ID fields available for an example."
            ),
            "",
            "## Expression-to-clinical sample matching",
            "",
            (
                "Sample IDs are compared exactly after trimming outer whitespace only. No barcode "
                "truncation, case conversion, or sample-type filtering is applied."
            ),
            "",
            "| Measure | Count |",
            "|---|---:|",
            f"| Expression sample columns | {expression.sample_columns:,} |",
            f"| Unique nonblank expression sample IDs | {matching.expression_unique_ids:,} |",
            f"| Clinical sample rows | {clinical.row_count:,} |",
            f"| Unique nonblank clinical sample IDs | {matching.clinical_unique_ids:,} |",
            f"| Duplicate exact clinical sample IDs | {len(clinical.duplicate_sample_ids):,} |",
            f"| Clinical rows with malformed width | {clinical.malformed_rows:,} |",
            f"| Exact sample-ID intersection | {intersection_size:,} |",
            f"| Expression-only unique IDs | {len(matching.expression_only_ids):,} |",
            f"| Clinical-only unique IDs | {len(matching.clinical_only_ids):,} |",
            (
                "| Expression IDs matched (intersection / unique expression IDs) | "
                f"{_format_percentage(intersection_size, expression_set_size)} |"
            ),
            (
                "| Clinical IDs matched (intersection / unique clinical IDs) | "
                f"{_format_percentage(intersection_size, clinical_set_size)} |"
            ),
            "",
            "Unmatched IDs are reported, not discarded:",
            "",
            f"- Expression-only: {_code_list(matching.expression_only_ids)}.",
            f"- Clinical-only: {_code_list(matching.clinical_only_ids)}.",
            f"- Clinical rows with missing sample IDs: **{clinical.missing_sample_ids}**.",
        ]
    )

    patient_count = len(patients.patient_sample_counts)
    distribution_text = (
        "; ".join(
            f"{sample_count} sample(s): {patient_count_for_count:,} patients"
            for sample_count, patient_count_for_count in (
                patients.samples_per_patient_distribution.items()
            )
        )
        or "No parseable TCGA patient IDs."
    )
    lines.extend(
        [
            "",
            "## Patient/sample relationship in expression",
            "",
            f"- Unique patients derived from expression sample barcodes: **{patient_count:,}**.",
            "- Patients with more than one expression sample: "
            f"**{patients.patients_with_multiple_samples:,}**.",
            f"- Distribution of expression samples per patient: {distribution_text}.",
            (
                "- Expression sample IDs that did not match the documented TCGA barcode "
                f"pattern: **{len(patients.unparseable_sample_ids)}**."
            ),
        ]
    )
    if patients.unparseable_sample_ids:
        lines.append(f"  - Examples: {_code_list(patients.unparseable_sample_ids[:10])}.")
    if patients.duplicate_patient_examples:
        lines.append(
            "- Examples of patients with multiple expression samples "
            "(no sample retained or discarded):"
        )
        for patient_id, sample_ids in patients.duplicate_patient_examples:
            lines.append(f"  - `{patient_id}`: {_code_list(sample_ids)}")
    else:
        lines.append("- No patients with more than one expression sample were detected.")

    lines.extend(
        [
            "",
            "## Clinical fields and missingness",
            "",
            (
                "Exact field names are discovered from the downloaded clinical header; categories "
                "are based on column-name patterns and are not a clinical data dictionary. "
                "Multiple receptor fields can reflect different source forms, assays, or "
                "specimen contexts; "
                "this audit does not select a preferred one."
            ),
            "",
        ]
    )
    for category, columns in clinical.key_columns.items():
        lines.append(f"### {category}")
        lines.append("")
        if columns:
            lines.append(_code_list(columns))
        else:
            lines.append("No matching column names are present in this file.")
        lines.append("")

    lines.extend(_render_missingness(clinical))
    lines.extend(
        [
            "",
            "## PAM50Call_RNAseq",
            "",
        ]
    )
    if clinical.pam50_column:
        labelled = sum(clinical.pam50_counts.values())
        lines.append(f"- Target field: `{clinical.pam50_column}`.")
        lines.append(
            f"- Nonmissing labels: **{labelled:,}**; missing labels: "
            f"**{clinical.pam50_missing:,} / {clinical.row_count:,}**."
        )
        lines.append("- Raw label counts (values preserved exactly; no normalization applied):")
        for label, count in sorted(clinical.pam50_counts.items()):
            lines.append(f"  - `{label}`: {count:,}")
        if clinical.pam50_case_or_spacing_variants:
            lines.append(
                "- Values differing only by case/outer whitespace (review before harmonizing):"
            )
            for group in clinical.pam50_case_or_spacing_variants:
                lines.append(f"  - {_code_list(group)}")
        else:
            lines.append(
                "- No within-field label pairs differing only by case or whitespace were detected."
            )
        if "PAM50_mRNA_nature2012" in clinical.pam50_other_field_counts:
            other = clinical.pam50_other_field_counts["PAM50_mRNA_nature2012"]
            other_values = {key: value for key, value in other.items() if key != "<MISSING>"}
            other_label_text = ", ".join(
                f"`{label}`: {count:,}" for label, count in sorted(other_values.items())
            )
            lines.append(
                "- A separate source field, `PAM50_mRNA_nature2012`, uses a different call "
                f"({sum(other_values.values()):,} nonmissing; "
                f"{other.get('<MISSING>', 0):,} missing): {other_label_text}. "
                "Do not merge it with `PAM50Call_RNAseq` without a documented decision."
            )
        lines.append(
            "- The selected field uses abbreviated source labels (for example, "
            "`LumA`, `Her2`, `Basal`, `Normal`). Possible display-name correspondences "
            "to review against project terminology are `LumA`/Luminal A, `LumB`/Luminal B, "
            "`Her2`/HER2-enriched, `Basal`/Basal-like, and `Normal`/Normal-like. "
            "These are not recodes and have not been applied."
        )
    else:
        lines.append(
            "The exact field `PAM50Call_RNAseq` was not found; no PAM50 distribution "
            "could be calculated."
        )

    lines.extend(
        [
            "",
            "## Important observations and limitations",
            "",
            (
                "- Expression and clinical sample counts differ; the exact intersection and all "
                "unmatched sample IDs are shown above. Matching is identifier-only, not evidence "
                "that labels are valid for every expression profile."
            ),
            (
                "- Patient-level grouping uses the first three TCGA barcode fields. The clinical "
                "field named `patient_id` is a short suffix where checked; do not treat it as a "
                "unique full barcode without confirmation."
            ),
            (
                "- The source contains multiple PAM50-related fields with different completeness "
                "and label vocabularies. `PAM50Call_RNAseq` is reported as-is; no label mapping or "
                "sample inclusion decision has been made."
            ),
            (
                "- The clinical file has no column whose name indicates race or ethnicity; "
                "demographic subgroup analysis needs a separately documented source/decision."
            ),
            (
                "- HiSeqV2 values are `log2(normalized_count + 1)`, not raw integer counts. They "
                "may be used for visualization or documented ML preprocessing, but must not be "
                "passed directly to DESeq2/PyDESeq2 as raw counts."
            ),
            (
                "- This is an identifier and metadata audit only. No genes or samples were "
                "filtered, no duplicate aliquot was selected, and no train/test split, "
                "differential expression, GSEA, or machine learning was performed."
            ),
            "",
            "## Decisions reserved for human review",
            "",
            (
                "1. Which PAM50 field and label vocabulary should define the later endpoint? This "
                "audit does not substitute `PAM50_mRNA_nature2012` for `PAM50Call_RNAseq`."
            ),
            (
                "2. Which ER/PR/HER2 fields should be authoritative if receptor-status analysis is "
                "later needed? Candidate fields have different missingness and contexts."
            ),
            (
                "3. How should unmatched samples and multiple samples per patient be handled "
                "later? No unmatched IDs or duplicate-patient samples were dropped here."
            ),
            (
                "4. A raw-count source is required before count-based differential expression; "
                "HiSeqV2 must remain labeled as normalized expression."
            ),
            "",
            "## Reproducibility",
            "",
            (
                "Regenerate this report from the repository root with "
                "`python scripts/audit_data.py`. The script reads the matrix one row at a time, "
                "creates no split, and writes this "
                "Markdown report; downloaded source files remain under Git-ignored `data/raw/`."
            ),
            "",
        ]
    )
    return "\n".join(lines)


def _audit_date() -> str:
    """Use Africa/Lagos date when available; otherwise use the machine's local date."""
    try:
        return datetime.now(ZoneInfo("Africa/Lagos")).date().isoformat()
    except ZoneInfoNotFoundError:
        return datetime.now().astimezone().date().isoformat()


def run_audit(
    expression_path: Path = EXPRESSION_PATH,
    clinical_path: Path = CLINICAL_PATH,
    report_path: Path | None = REPORT_PATH,
    audit_date: str | None = None,
) -> AuditResult:
    """Audit both files and optionally write the generated Markdown report."""
    expression = inspect_expression(expression_path)
    clinical = inspect_clinical(clinical_path)
    matching = match_samples(expression, clinical)
    patients = audit_expression_patients(expression.sample_ids)
    markdown = _render_report(
        expression=expression,
        clinical=clinical,
        matching=matching,
        patients=patients,
        audit_date=audit_date or _audit_date(),
    )

    if report_path is not None:
        report_path = Path(report_path)
        report_path.parent.mkdir(parents=True, exist_ok=True)
        report_path.write_text(markdown, encoding="utf-8")

    return AuditResult(
        expression=expression,
        clinical=clinical,
        matching=matching,
        patients=patients,
        markdown=markdown,
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Audit TCGA-BRCA expression and clinical identifiers without preprocessing."
    )
    parser.add_argument(
        "--expression",
        type=Path,
        default=EXPRESSION_PATH,
        help="Path to the downloaded HiSeqV2 gzip matrix.",
    )
    parser.add_argument(
        "--clinical",
        type=Path,
        default=CLINICAL_PATH,
        help="Path to the downloaded BRCA clinical matrix.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=REPORT_PATH,
        help=f"Markdown report destination (default: {REPORT_PATH.relative_to(PROJECT_ROOT)}).",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        result = run_audit(
            expression_path=args.expression,
            clinical_path=args.clinical,
            report_path=args.output,
        )
    except (OSError, EOFError, csv.Error, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1

    print(f"Expression samples: {result.expression.sample_columns:,}")
    print(f"Expression genes: {result.expression.gene_rows:,}")
    print(f"Clinical rows: {result.clinical.row_count:,}")
    print(f"Clinical columns: {len(result.clinical.columns):,}")
    print(f"Matched unique sample IDs: {len(result.matching.intersection_ids):,}")
    print(f"Report written: {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
