"""Build a metadata-only one-sample-per-patient TCGA-BRCA PAM50 cohort.

Run from the repository root with:

    python scripts/build_cohort.py

Only the expression matrix header is read. Expression values are never loaded.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import sys
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path
from zoneinfo import ZoneInfo, ZoneInfoNotFoundError

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_PATH = PROJECT_ROOT / "src"
if str(SRC_PATH) not in sys.path:
    sys.path.insert(0, str(SRC_PATH))

from oncorna.cohort import (  # noqa: E402
    PAM50_CLASS_ORDER,
    CohortBuildResult,
    SampleAuditRecord,
    build_classification_cohort,
)

EXPRESSION_PATH = PROJECT_ROOT / "data" / "raw" / "TCGA-BRCA_HiSeqV2.tsv.gz"
CLINICAL_PATH = PROJECT_ROOT / "data" / "raw" / "TCGA-BRCA_BRCA_clinicalMatrix.tsv"
OUTPUT_PATH = PROJECT_ROOT / "data" / "processed" / "classification_cohort.tsv"
AUDIT_PATH = PROJECT_ROOT / "data" / "processed" / "cohort_sample_audit.tsv"
REPORT_PATH = PROJECT_ROOT / "docs" / "cohort_definition.md"

REQUIRED_CLINICAL_FIELDS = {
    "sampleID",
    "PAM50Call_RNAseq",
    "sample_type",
    "sample_type_id",
    "_PATIENT",
}
OPTIONAL_CLINICAL_FIELDS = {
    "patient_id",
    "bcr_patient_barcode",
    "bcr_sample_barcode",
}
COHORT_FIELDS = (
    "sample_id",
    "patient_id",
    "pam50_original_label",
    "pam50_normalized_label",
    "sample_type",
    "sample_type_id",
    "sample_type_code_from_barcode",
    "sample_retention_reason",
)
AUDIT_FIELDS = (
    "sample_id",
    "patient_id",
    "clinical_match_status",
    "pam50_original_label",
    "pam50_normalized_label",
    "pam50_label_status",
    "sample_type",
    "sample_type_id",
    "sample_type_code_from_barcode",
    "sample_type_classification",
    "sample_type_evidence_status",
    "patient_barcode_validation_status",
    "clinical_patient_id_status",
    "clinical_bcr_patient_barcode_status",
    "bcr_sample_barcode_status",
    "selected_sample_id_for_patient",
    "is_selected_patient_sample",
    "is_in_classification_cohort",
    "sample_retention_reason",
)
DISPOSITION_DESCRIPTIONS = {
    "retained_primary_tumor": "Retained: selected primary tumor with a mapped PAM50 label",
    "retained_metastatic": "Retained: selected metastatic sample with a mapped PAM50 label",
    "retained_solid_tissue_normal": (
        "Retained: selected solid-tissue-normal sample with a mapped PAM50 label"
    ),
    "excluded_missing_pam50": "Excluded: selected preferred sample has missing PAM50",
    "unmapped_pam50_label_requires_review": (
        "Review: selected sample has a nonmissing, unmapped PAM50 label"
    ),
    "not_selected_lower_priority_sample_type": (
        "Not selected: a higher-priority sample type exists for this patient"
    ),
    "not_selected_lexicographic_tie_break": (
        "Not selected: lexicographically later ID among same-priority samples"
    ),
    "no_matching_clinical_metadata": "Excluded: no exact clinical sampleID match",
    "duplicate_clinical_sample_id_requires_review": (
        "Review: clinical sampleID occurs more than once"
    ),
    "unparseable_tcga_patient_id": "Review: patient ID could not be derived from TCGA barcode",
    "patient_identifier_conflict_requires_review": (
        "Review: clinical _PATIENT conflicts with sample barcode"
    ),
    "sample_type_conflict_requires_review": (
        "Review: sample-type evidence is conflicting or unmapped"
    ),
}
SAMPLE_TYPE_DISPLAY_ORDER = ("01", "06", "11")


def read_expression_sample_ids(path: Path) -> list[str]:
    """Read only the expression header; never load expression values."""
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"Expression matrix not found: {path}")
    try:
        with gzip.open(path, "rt", encoding="utf-8-sig", newline="") as handle:
            header = next(csv.reader(handle, delimiter="\t"), None)
    except (OSError, EOFError, csv.Error) as error:
        raise ValueError(f"Could not read expression header from {path}: {error}") from error
    if not header or len(header) < 2:
        raise ValueError(f"Expression matrix has no usable sample header: {path}")
    sample_ids = [sample.strip() for sample in header[1:]]
    if any(not sample_id for sample_id in sample_ids):
        raise ValueError("Expression header contains a blank sample ID")
    duplicate_ids = [sample_id for sample_id, count in Counter(sample_ids).items() if count > 1]
    if duplicate_ids:
        examples = ", ".join(sorted(duplicate_ids)[:5])
        raise ValueError(f"Expression header contains duplicate sample IDs: {examples}")
    return sample_ids


def read_clinical_metadata(path: Path) -> list[dict[str, str | None]]:
    """Read only the clinical fields needed for sample, label, and ID checks."""
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"Clinical matrix not found: {path}")
    selected_fields = REQUIRED_CLINICAL_FIELDS | OPTIONAL_CLINICAL_FIELDS
    try:
        with path.open("rt", encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            columns = reader.fieldnames
            if not columns:
                raise ValueError(f"Clinical matrix has no usable header: {path}")
            missing_fields = REQUIRED_CLINICAL_FIELDS - set(columns)
            if missing_fields:
                missing = ", ".join(sorted(missing_fields))
                raise ValueError(f"Clinical matrix is missing required fields: {missing}")
            if len(columns) != len(set(columns)):
                raise ValueError("Clinical matrix contains duplicate column names")

            rows: list[dict[str, str | None]] = []
            for row_number, row in enumerate(reader, start=2):
                if None in row or any(row.get(column) is None for column in columns):
                    raise ValueError(f"Clinical matrix row {row_number} has a width mismatch")
                rows.append(
                    {field: row.get(field) for field in selected_fields if field in columns}
                )
    except (OSError, csv.Error) as error:
        raise ValueError(f"Could not read clinical metadata from {path}: {error}") from error
    return rows


def _write_tsv(path: Path, fieldnames: tuple[str, ...], rows: list[dict[str, str]]) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _record_to_cohort_row(record: SampleAuditRecord) -> dict[str, str]:
    return {
        "sample_id": record.sample_id,
        "patient_id": record.patient_id or "",
        "pam50_original_label": record.pam50_original_label,
        "pam50_normalized_label": record.pam50_normalized_label or "",
        "sample_type": record.sample_type,
        "sample_type_id": record.sample_type_id,
        "sample_type_code_from_barcode": record.sample_type_code_from_barcode or "",
        "sample_retention_reason": record.sample_retention_reason,
    }


def _record_to_audit_row(record: SampleAuditRecord) -> dict[str, str]:
    row = _record_to_cohort_row(record)
    row.update(
        {
            "clinical_match_status": record.clinical_match_status,
            "pam50_label_status": record.pam50_label_status,
            "sample_type_classification": record.sample_type_classification,
            "sample_type_evidence_status": record.sample_type_evidence_status,
            "patient_barcode_validation_status": record.patient_barcode_validation_status,
            "clinical_patient_id_status": record.clinical_patient_id_status,
            "clinical_bcr_patient_barcode_status": record.clinical_bcr_patient_barcode_status,
            "bcr_sample_barcode_status": record.bcr_sample_barcode_status,
            "selected_sample_id_for_patient": record.selected_sample_id_for_patient,
            "is_selected_patient_sample": str(record.is_selected_patient_sample).lower(),
            "is_in_classification_cohort": str(record.is_in_classification_cohort).lower(),
        }
    )
    return row


def _today() -> str:
    try:
        return datetime.now(ZoneInfo("Africa/Lagos")).date().isoformat()
    except ZoneInfoNotFoundError:
        return datetime.now().astimezone().date().isoformat()


def _counter_for(records: list[SampleAuditRecord], field: str) -> Counter[str]:
    return Counter(str(getattr(record, field)) for record in records)


def _identifier_check_row(
    label: str,
    counts: Counter[str],
    *,
    match_key: str = "match",
) -> str:
    return (
        f"| {label} | {counts.get(match_key, 0):,} | "
        f"{counts.get('mismatch', 0):,} | {counts.get('missing', 0):,} |"
    )


def _render_cohort_definition(result: CohortBuildResult, generated_date: str) -> str:
    matched = [
        record for record in result.sample_records if record.clinical_match_status == "matched"
    ]
    type_counts = Counter(record.sample_type_code_from_barcode or "unparsed" for record in matched)
    type_data: dict[str, dict[str, object]] = defaultdict(dict)
    for code in type_counts:
        records = [record for record in matched if record.sample_type_code_from_barcode == code]
        type_data[code] = {
            "label": records[0].sample_type if records else "Unknown",
            "count": len(records),
            "missing": sum(record.pam50_label_status == "missing" for record in records),
        }

    patient_type_groups: dict[str, list[str]] = defaultdict(list)
    for record in matched:
        if record.patient_id:
            patient_type_groups[record.patient_id].append(
                record.sample_type_code_from_barcode or "unparsed"
            )
    composition_counts = Counter(
        "+".join(sorted(Counter(codes).elements()))
        for codes in patient_type_groups.values()
        if len(codes) > 1
    )

    normalized_counts = Counter(record.pam50_normalized_label for record in result.cohort)
    expression_patient_count = sum(result.patient_sample_count_distribution.values())
    patient_sample_distribution = _distribution_text(result.patient_sample_count_distribution)
    secondary_missing = sum(
        record.pam50_label_status == "missing"
        and record.sample_retention_reason == "not_selected_lower_priority_sample_type"
        for record in result.sample_records
    )
    primary_selected_count = sum(
        record.is_selected_patient_sample and record.sample_type_category == "primary_tumor"
        for record in result.sample_records
    )
    patient_id_checks = _counter_for(matched, "patient_barcode_validation_status")
    short_patient_checks = _counter_for(matched, "clinical_patient_id_status")
    full_patient_checks = _counter_for(matched, "clinical_bcr_patient_barcode_status")
    sample_barcode_checks = _counter_for(matched, "bcr_sample_barcode_status")
    sample_type_evidence = _counter_for(matched, "sample_type_evidence_status")
    sample_type_evidence_text = ", ".join(
        f"`{key}` {value:,}" for key, value in sorted(sample_type_evidence.items())
    )

    ambiguous_ids = sorted(
        {
            record.sample_id
            for record in matched
            if record.clinical_patient_id_status == "mismatch"
            or record.clinical_bcr_patient_barcode_status == "mismatch"
            or record.bcr_sample_barcode_status == "mismatch"
        }
    )
    ambiguous_id_text = ", ".join(f"`{sample_id}`" for sample_id in ambiguous_ids) or "None"

    lines = [
        "# TCGA-BRCA classification cohort definition",
        "",
        f"**Generated:** {generated_date}<br>",
        "**Status:** Phase 1C cohort definition only. No expression QC, gene filtering, splitting, "
        "differential expression, GSEA, or machine learning was performed.",
        "",
        "## Cohort rules",
        "",
        "The primary outcome is the clinical field `PAM50Call_RNAseq`. The source value is copied "
        "unchanged to `pam50_original_label`; the separate normalized field maps `LumA` → "
        "`Luminal A`, `LumB` → `Luminal B`, `Basal` → `Basal-like`, `Her2` → "
        "`HER2-enriched`, and `Normal` → `Normal-like`. All five classes are retained.",
        "",
        "Sample type is read from the fourth TCGA barcode field (the first two digits are the "
        "sample-type code) and checked against the clinical `sample_type_id` and `sample_type` "
        "fields. GDC defines code 01 as Primary Solid Tumor, 06 as Metastatic, and 11 as Solid "
        "Tissue Normal ([sample-type code table](https://gdc.cancer.gov/resources-tcga-users/"
        "tcga-code-tables/sample-type-codes); [TCGA barcode guide](https://docs.gdc.cancer.gov/"
        "Encyclopedia/pages/TCGA_Barcode/)). The downloaded clinical text labels for those codes "
        "are `Primary Tumor`, `Metastatic`, and `Solid Tissue Normal`.",
        "",
        "The matrix also contains `bcr_sample_barcode`, `tumor_tissue_site`, "
        "`tissue_source_site`, `tissue_prospective_collection_indicator`, and "
        "`tissue_retrospective_collection_indicator`. These describe barcode, anatomic site, "
        "source site, or collection context; they do not replace the explicit sample-type "
        "code/text for distinguishing primary tumor, metastasis, and solid-tissue normal samples.",
        "",
        "For each participant derived from the first three sample-barcode fields, the exact "
        "selection hierarchy is: **01 Primary Tumor → 06 Metastatic → 11 Solid Tissue Normal**. "
        "Within a type-priority tie, the lexicographically smallest full sample ID wins. No "
        "expression values, PAM50 labels/confidence, survival, or downstream results are used to "
        "choose a sample.",
        "",
        "One sample is selected per patient before applying PAM50 availability. If the selected "
        "sample has missing PAM50, that patient is excluded; the builder does not fall back to a "
        "secondary sample with a nonmissing label. This preserves the primary-tumor preference. "
        "Unrecognized or conflicting sample-type evidence and conflicting `_PATIENT` values are "
        "flagged for human review rather than ranked automatically.",
        "",
        "## Starting data and identifier checks",
        "",
        f"- Expression sample columns at start: **{result.expression_sample_count:,}**.",
        (
            f"- Clinical sample rows: **{result.clinical_row_count:,}**; "
            f"unique nonblank clinical `sampleID`s: **{result.clinical_unique_sample_count:,}**."
        ),
        (
            f"- Exact expression/clinical `sampleID` matches: "
            f"**{result.matched_expression_count:,} / {result.expression_sample_count:,}**."
        ),
        (
            "- Expression samples without an exact clinical match: "
            f"**{result.expression_only_sample_count:,}**."
        ),
        (
            f"- Clinical-only sample IDs (not in expression): "
            f"**{result.clinical_only_sample_count:,}**; these are outside the "
            "expression-sample denominator."
        ),
        f"- Duplicate clinical sample IDs: **{len(result.duplicate_clinical_sample_ids):,}**.",
        (
            "- Unique expression participants from TCGA barcode parsing: "
            f"**{expression_patient_count:,}**."
        ),
        (
            "- Patients with multiple expression samples: "
            f"**{result.multiple_sample_patient_count:,}**."
        ),
        "",
        "### Expression sample types",
        "",
        (
            "| Barcode code | Clinical `sample_type` | Matched expression samples | "
            "Missing PAM50 among these samples |"
        ),
        "|---|---|---:|---:|",
    ]
    for code in (
        *SAMPLE_TYPE_DISPLAY_ORDER,
        *sorted(set(type_data) - set(SAMPLE_TYPE_DISPLAY_ORDER)),
    ):
        if code not in type_data:
            continue
        item = type_data[code]
        lines.append(f"| `{code}` | `{item['label']}` | {item['count']:,} | {item['missing']:,} |")

    lines.extend(
        [
            "",
            "For all matched expression IDs, the barcode code, clinical `sample_type_id`, and "
            "clinical `sample_type` text agree. Every derived participant ID also agrees with "
            "clinical `_PATIENT` in the matched rows. The older `patient_id` and "
            "`bcr_patient_barcode` fields are not used as the grouping key because they contain "
            "missing and truncated values (details below).",
            "",
            "### TCGA identifier cross-checks on matched expression rows",
            "",
            "| Clinical field/check | Matches | Mismatches | Missing |",
            "|---|---:|---:|---:|",
            _identifier_check_row("Derived patient ID vs `_PATIENT`", patient_id_checks),
            _identifier_check_row("Derived patient suffix vs `patient_id`", short_patient_checks),
            _identifier_check_row(
                "Derived patient ID vs `bcr_patient_barcode`", full_patient_checks
            ),
            _identifier_check_row(
                "Derived patient and sample type vs `bcr_sample_barcode`",
                sample_barcode_checks,
                match_key="patient_and_type_match",
            ),
            "",
            "Two expression rows (`TCGA-E9-A1NA-01` and `TCGA-E9-A1NA-11`) for participant "
            "`TCGA-E9-A1NA` have `patient_id`=`A1` and `bcr_patient_barcode`=`TCGA-E9-A1`, "
            "while `sampleID` and `_PATIENT` identify "
            "`TCGA-E9-A1NA`. The same two `bcr_sample_barcode` values use the truncated "
            "participant portion. `_PATIENT`, barcode sample-type code, `sample_type_id`, and "
            "`sample_type` agree for both rows. These rows are grouped using the sample ID's "
            "first three fields, cross-checked with `_PATIENT`; the truncated legacy fields are "
            "not used for grouping. Human review of this source inconsistency is still warranted.",
            "",
            f"- Sample IDs with conflicting legacy barcode fields: {ambiguous_id_text}.",
            "",
            f"Sample-type evidence status across matched expression rows: "
            f"{sample_type_evidence_text}.",
            "",
            "## One-patient-one-sample resolution",
            "",
            (
                f"- The expression sample set represents **{expression_patient_count:,}** "
                f"participants; **{result.multiple_sample_patient_count:,}** have multiple samples."
            ),
            f"- Samples per patient before selection: {patient_sample_distribution}.",
            (
                "- Participants with exactly one observed type-01 primary tumor sample: "
                f"**{primary_selected_count:,}**."
            ),
            "- Type compositions among multi-sample participants:",
        ]
    )
    if composition_counts:
        for composition, count in sorted(composition_counts.items()):
            lines.append(f"  - `{composition}`: {count:,} participant(s)")
    else:
        lines.append("  - No multi-sample participants.")

    lines.extend(
        [
            "",
            "The 121 non-primary samples are not selected because every represented participant "
            "has a type-01 primary tumor sample. Of these 121 secondary samples, 112 have a "
            "nonmissing PAM50 label and 9 are missing it; neither label availability nor label "
            "identity changes the primary-tumor selection.",
            "",
            "## Exclusions and final cohort",
            "",
            (
                f"Among all **{result.matched_expression_count:,}** expression samples with "
                f"matching clinical metadata, **{result.missing_pam50_matched_sample_count:,}** "
                "have missing `PAM50Call_RNAseq` across all sample types. After the sample-type "
                "hierarchy selects one sample per patient, "
                f"**{result.missing_pam50_selected_sample_count:,}** selected primary-tumor "
                "samples are excluded for missing PAM50. The other "
                f"**{secondary_missing:,}** missing labels are on secondary samples already "
                "superseded by the primary-tumor rule and are not double-counted."
            ),
            "",
            "Expression-sample disposition (counts sum to the starting expression sample count):",
            "",
            "| Reason | Count |",
            "|---|---:|",
        ]
    )
    disposition_counts = result.disposition_counts
    ordered_reasons = list(DISPOSITION_DESCRIPTIONS)
    ordered_reasons.extend(sorted(set(disposition_counts) - set(ordered_reasons)))
    for reason in ordered_reasons:
        description = DISPOSITION_DESCRIPTIONS.get(reason, reason.replace("_", " "))
        lines.append(f"| {description} (`{reason}`) | {disposition_counts.get(reason, 0):,} |")
    lines.append("")
    lines.append(
        f"- Clinical-only sample IDs outside the expression denominator: "
        f"**{result.clinical_only_sample_count:,}**."
    )
    lines.extend(
        [
            "",
            "- Selected one-sample-per-patient records before PAM50 eligibility: "
            f"**{result.selected_patient_sample_count:,}**.",
            f"- Final classification cohort: **{len(result.cohort):,} unique patients**.",
            "",
            "### Final normalized PAM50 distribution",
            "",
            "| Normalized PAM50 label | Samples/patients |",
            "|---|---:|",
        ]
    )
    for label in PAM50_CLASS_ORDER:
        lines.append(f"| `{label}` | {normalized_counts.get(label, 0):,} |")

    lines.extend(
        [
            "",
            "All five classes are retained, including `Normal-like`; no minimum class-size rule "
            "has been applied. `classification_cohort.tsv` contains metadata only. "
            "`cohort_sample_audit.tsv` contains one disposition row per expression sample so "
            "every selection and exclusion is traceable. Both files are written under the "
            "Git-ignored `data/processed/` directory; the full expression matrix is not copied.",
            "",
            "## Reproducibility and scope boundary",
            "",
            "Regenerate the cohort tables and this document from the repository root with "
            "`python scripts/build_cohort.py`. The builder reads the expression header only, "
            "matches clinical `sampleID`s exactly after outer-whitespace trimming, and writes "
            "metadata only. It does not perform expression QC, gene filtering, train/test "
            "splitting, differential expression, GSEA, or machine learning.",
            "",
            "## Human review items",
            "",
            "1. Review the truncated `TCGA-E9-A1NA` legacy clinical barcode fields described "
            "above; the grouping key is supported by `_PATIENT` and expression `sampleID`.",
            "2. Confirm that selecting the primary tumor before PAM50 filtering, with no fallback "
            "when that selected label is missing, remains the desired rule for later work.",
            "3. Confirm any future sample-type codes outside 01/06/11 before assigning them a "
            "selection rank; this build does not guess an order for unrecognized codes.",
            "",
        ]
    )
    return "\n".join(lines)


def _distribution_text(distribution: dict[int, int]) -> str:
    if not distribution:
        return "No parseable TCGA participant IDs."
    return "; ".join(
        f"{sample_count} sample(s): {patient_count:,} patient(s)"
        for sample_count, patient_count in sorted(distribution.items())
    )


def run_build(
    expression_path: Path = EXPRESSION_PATH,
    clinical_path: Path = CLINICAL_PATH,
    cohort_output: Path = OUTPUT_PATH,
    audit_output: Path = AUDIT_PATH,
    report_output: Path = REPORT_PATH,
) -> CohortBuildResult:
    """Read source identifiers/metadata and write cohort, audit, and report files."""
    expression_sample_ids = read_expression_sample_ids(expression_path)
    clinical_rows = read_clinical_metadata(clinical_path)
    result = build_classification_cohort(expression_sample_ids, clinical_rows)

    _write_tsv(
        cohort_output,
        COHORT_FIELDS,
        [_record_to_cohort_row(record) for record in result.cohort],
    )
    _write_tsv(
        audit_output,
        AUDIT_FIELDS,
        [_record_to_audit_row(record) for record in result.sample_records],
    )
    report_output = Path(report_output)
    report_output.parent.mkdir(parents=True, exist_ok=True)
    report_output.write_text(
        _render_cohort_definition(result, _today()),
        encoding="utf-8",
    )
    return result


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Build the TCGA-BRCA one-patient-one-sample PAM50 cohort metadata."
    )
    parser.add_argument("--expression", type=Path, default=EXPRESSION_PATH)
    parser.add_argument("--clinical", type=Path, default=CLINICAL_PATH)
    parser.add_argument("--cohort-output", type=Path, default=OUTPUT_PATH)
    parser.add_argument("--audit-output", type=Path, default=AUDIT_PATH)
    parser.add_argument("--report", type=Path, default=REPORT_PATH)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        result = run_build(
            expression_path=args.expression,
            clinical_path=args.clinical,
            cohort_output=args.cohort_output,
            audit_output=args.audit_output,
            report_output=args.report,
        )
    except (OSError, EOFError, csv.Error, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1

    print(f"Expression samples: {result.expression_sample_count:,}")
    print(f"Matched expression/clinical samples: {result.matched_expression_count:,}")
    print(f"Patients with multiple samples: {result.multiple_sample_patient_count:,}")
    print(f"Final unique-patient cohort: {len(result.cohort):,}")
    print(f"Cohort table: {args.cohort_output}")
    print(f"Sample audit: {args.audit_output}")
    print(f"Cohort definition: {args.report}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
