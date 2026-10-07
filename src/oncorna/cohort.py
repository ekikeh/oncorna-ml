"""Identifier-safe cohort construction for the TCGA-BRCA PAM50 task.

This module operates on sample identifiers and clinical metadata only. It never
loads or transforms expression values.
"""

from __future__ import annotations

import re
from collections import Counter, defaultdict
from dataclasses import dataclass
from typing import Iterable, Mapping, Sequence

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

PAM50_LABEL_MAP = {
    "LumA": "Luminal A",
    "LumB": "Luminal B",
    "Basal": "Basal-like",
    "Her2": "HER2-enriched",
    "Normal": "Normal-like",
}
PAM50_CLASS_ORDER = (
    "Luminal A",
    "Luminal B",
    "Basal-like",
    "HER2-enriched",
    "Normal-like",
)

TCGA_CENTER_CODE = re.compile(r"^[A-Z0-9]{2}$", re.IGNORECASE)
TCGA_PARTICIPANT_CODE = re.compile(r"^[A-Z0-9]{4}$", re.IGNORECASE)
TCGA_SAMPLE_SEGMENT = re.compile(r"^(?P<code>\d{2})(?:[A-Z])?$", re.IGNORECASE)


@dataclass(frozen=True)
class _SampleTypeDefinition:
    code: str
    label: str
    category: str
    priority: int
    aliases: tuple[str, ...]


_SAMPLE_TYPES = {
    "01": _SampleTypeDefinition(
        code="01",
        label="Primary Tumor",
        category="primary_tumor",
        priority=0,
        aliases=("primary tumor", "primary solid tumor"),
    ),
    "06": _SampleTypeDefinition(
        code="06",
        label="Metastatic",
        category="metastatic",
        priority=1,
        aliases=("metastatic",),
    ),
    "11": _SampleTypeDefinition(
        code="11",
        label="Solid Tissue Normal",
        category="solid_tissue_normal",
        priority=2,
        aliases=("solid tissue normal",),
    ),
}
_SAMPLE_TYPE_TEXT_TO_CODE = {
    alias: code for code, definition in _SAMPLE_TYPES.items() for alias in definition.aliases
}


@dataclass(frozen=True)
class SampleTypeClassification:
    """Sample-type interpretation plus the consistency of its source evidence."""

    code: str | None
    label: str
    category: str
    priority: int | None
    evidence_status: str

    @property
    def is_resolved(self) -> bool:
        return self.priority is not None and self.evidence_status in {
            "consistent",
            "partial_but_consistent",
        }


@dataclass
class SampleAuditRecord:
    """Metadata-only disposition for one expression sample column."""

    sample_id: str
    patient_id: str | None
    clinical_match_status: str
    pam50_original_label: str
    pam50_normalized_label: str | None
    pam50_label_status: str
    sample_type: str
    sample_type_id: str
    sample_type_code_from_barcode: str | None
    sample_type_classification: str
    sample_type_category: str
    sample_type_priority: int | None
    sample_type_evidence_status: str
    patient_barcode_validation_status: str
    clinical_patient_id_status: str
    clinical_bcr_patient_barcode_status: str
    bcr_sample_barcode_status: str
    selected_sample_id_for_patient: str = ""
    is_selected_patient_sample: bool = False
    is_in_classification_cohort: bool = False
    sample_retention_reason: str = ""


@dataclass
class CohortBuildResult:
    """Cohort rows and explicit sample-level dispositions."""

    sample_records: list[SampleAuditRecord]
    cohort: list[SampleAuditRecord]
    expression_sample_count: int
    clinical_row_count: int
    clinical_unique_sample_count: int
    clinical_only_sample_count: int
    duplicate_clinical_sample_ids: list[str]
    multiple_sample_patient_count: int
    patient_sample_count_distribution: dict[int, int]

    @property
    def matched_expression_count(self) -> int:
        return sum(record.clinical_match_status == "matched" for record in self.sample_records)

    @property
    def expression_only_sample_count(self) -> int:
        return sum(record.clinical_match_status == "no_match" for record in self.sample_records)

    @property
    def selected_patient_sample_count(self) -> int:
        return sum(record.is_selected_patient_sample for record in self.sample_records)

    @property
    def missing_pam50_matched_sample_count(self) -> int:
        return sum(
            record.clinical_match_status == "matched" and record.pam50_label_status == "missing"
            for record in self.sample_records
        )

    @property
    def missing_pam50_selected_sample_count(self) -> int:
        return sum(
            record.is_selected_patient_sample and record.pam50_label_status == "missing"
            for record in self.sample_records
        )

    @property
    def disposition_counts(self) -> dict[str, int]:
        return dict(
            sorted(
                Counter(record.sample_retention_reason for record in self.sample_records).items()
            )
        )


def _is_missing(value: str | None) -> bool:
    return value is None or value.strip().casefold() in MISSING_VALUES


def derive_patient_id(sample_id: str) -> str | None:
    """Derive a TCGA participant ID from the first three barcode fields.

    ``TCGA-AR-A5QQ-01`` becomes ``TCGA-AR-A5QQ``. The source sample ID is not
    modified; its sample-type, vial, or aliquot suffix is omitted only from the
    derived grouping key.
    """
    if not sample_id:
        return None
    parts = sample_id.strip().split("-")
    if len(parts) < 3 or parts[0].upper() != "TCGA":
        return None
    if not TCGA_CENTER_CODE.fullmatch(parts[1]):
        return None
    if not TCGA_PARTICIPANT_CODE.fullmatch(parts[2]):
        return None
    return "-".join(parts[:3])


def parse_sample_type_code(sample_id: str) -> str | None:
    """Extract the two-digit TCGA sample-type code from barcode field four.

    The optional trailing letter is the vial identifier (for example, ``01A``);
    it is not part of the sample-type code.
    """
    if not sample_id:
        return None
    parts = sample_id.strip().split("-")
    if len(parts) < 4 or parts[0].upper() != "TCGA":
        return None
    match = TCGA_SAMPLE_SEGMENT.fullmatch(parts[3])
    return match.group("code") if match else None


def classify_sample_type(
    sample_id: str,
    sample_type: str | None = None,
    sample_type_id: str | None = None,
) -> SampleTypeClassification:
    """Classify a sample from barcode structure and available clinical fields.

    The current project mapping covers observed codes 01, 06, and 11. All
    non-missing evidence sources must agree. Unknown codes or inconsistent
    metadata receive no ranking and require review rather than a guessed type.
    """
    barcode_code = parse_sample_type_code(sample_id)
    evidence: list[tuple[str, str]] = []
    if barcode_code:
        evidence.append(("barcode", barcode_code))

    raw_type_id = (sample_type_id or "").strip()
    has_type_id = not _is_missing(raw_type_id)
    if has_type_id:
        if not re.fullmatch(r"\d{2}", raw_type_id):
            return SampleTypeClassification(
                code=barcode_code,
                label="Conflicting sample type metadata",
                category="ambiguous",
                priority=None,
                evidence_status="invalid_sample_type_id",
            )
        evidence.append(("sample_type_id", raw_type_id))

    normalized_text = " ".join((sample_type or "").strip().casefold().split())
    has_type_text = not _is_missing(normalized_text)
    if has_type_text:
        text_code = _SAMPLE_TYPE_TEXT_TO_CODE.get(normalized_text)
        if text_code is None:
            return SampleTypeClassification(
                code=barcode_code,
                label="Conflicting sample type metadata",
                category="ambiguous",
                priority=None,
                evidence_status="unrecognized_sample_type_label",
            )
        evidence.append(("sample_type", text_code))

    if not evidence:
        return SampleTypeClassification(
            code=None,
            label="Unknown sample type",
            category="unknown",
            priority=None,
            evidence_status="no_sample_type_evidence",
        )

    codes = {code for _, code in evidence}
    if len(codes) != 1:
        return SampleTypeClassification(
            code=barcode_code,
            label="Conflicting sample type metadata",
            category="ambiguous",
            priority=None,
            evidence_status="sample_type_sources_disagree",
        )

    code = codes.pop()
    definition = _SAMPLE_TYPES.get(code)
    if definition is None:
        return SampleTypeClassification(
            code=code,
            label=f"Unclassified sample type ({code})",
            category="unknown",
            priority=None,
            evidence_status="unmapped_sample_type_code",
        )

    all_sources_present = bool(barcode_code and has_type_id and has_type_text)
    return SampleTypeClassification(
        code=code,
        label=definition.label,
        category=definition.category,
        priority=definition.priority,
        evidence_status="consistent" if all_sources_present else "partial_but_consistent",
    )


def normalize_pam50(label: str | None) -> str | None:
    """Map the five observed source labels without modifying the original value.

    Outer whitespace is ignored for lookup, but case and spelling are not
    guessed. Missing values and non-missing labels outside the explicit map
    return ``None``; callers can distinguish them with ``_is_missing``.
    """
    if _is_missing(label):
        return None
    return PAM50_LABEL_MAP.get(label.strip())


def _compare_identifier(value: str | None, expected: str) -> str:
    if _is_missing(value):
        return "missing"
    return "match" if value.strip() == expected else "mismatch"


def _bcr_sample_barcode_status(
    barcode: str | None,
    patient_id: str | None,
    sample_type_code: str | None,
) -> str:
    if _is_missing(barcode):
        return "missing"
    bcr_patient_id = derive_patient_id(barcode or "")
    bcr_sample_type = parse_sample_type_code(barcode or "")
    if bcr_patient_id == patient_id and bcr_sample_type == sample_type_code:
        return "patient_and_type_match"
    return "mismatch"


def _make_record(
    sample_id: str,
    patient_id: str | None,
    row: Mapping[str, str | None] | None,
    clinical_match_status: str,
) -> SampleAuditRecord:
    row = row or {}
    original_label = row.get("PAM50Call_RNAseq") or ""
    normalized_label = normalize_pam50(original_label)
    if _is_missing(original_label):
        pam50_status = "missing"
    elif normalized_label is None:
        pam50_status = "unmapped"
    else:
        pam50_status = "mapped"

    source_sample_type = row.get("sample_type") or ""
    source_sample_type_id = row.get("sample_type_id") or ""
    type_result = classify_sample_type(
        sample_id,
        sample_type=source_sample_type,
        sample_type_id=source_sample_type_id,
    )
    barcode_code = parse_sample_type_code(sample_id)

    if row:
        full_patient_barcode = row.get("_PATIENT")
        patient_validation = (
            "not_comparable"
            if patient_id is None
            else _compare_identifier(full_patient_barcode, patient_id)
        )
        clinical_patient_id = row.get("patient_id")
        expected_short_id = patient_id.rsplit("-", maxsplit=1)[-1] if patient_id else ""
        patient_id_status = (
            "not_comparable"
            if patient_id is None
            else _compare_identifier(clinical_patient_id, expected_short_id)
        )
        bcr_patient_barcode_status = (
            "not_comparable"
            if patient_id is None
            else _compare_identifier(row.get("bcr_patient_barcode"), patient_id)
        )
        bcr_sample_status = _bcr_sample_barcode_status(
            row.get("bcr_sample_barcode"), patient_id, barcode_code
        )
    else:
        patient_validation = "not_checked"
        patient_id_status = "not_checked"
        bcr_patient_barcode_status = "not_checked"
        bcr_sample_status = "not_checked"

    return SampleAuditRecord(
        sample_id=sample_id,
        patient_id=patient_id,
        clinical_match_status=clinical_match_status,
        pam50_original_label=original_label,
        pam50_normalized_label=normalized_label,
        pam50_label_status=pam50_status,
        sample_type=source_sample_type,
        sample_type_id=source_sample_type_id,
        sample_type_code_from_barcode=barcode_code,
        sample_type_classification=type_result.label,
        sample_type_category=type_result.category,
        sample_type_priority=type_result.priority,
        sample_type_evidence_status=type_result.evidence_status,
        patient_barcode_validation_status=patient_validation,
        clinical_patient_id_status=patient_id_status,
        clinical_bcr_patient_barcode_status=bcr_patient_barcode_status,
        bcr_sample_barcode_status=bcr_sample_status,
    )


def select_one_sample_per_patient(
    samples: Sequence[SampleAuditRecord],
) -> dict[str, SampleAuditRecord]:
    """Choose one metadata-matched sample per patient by type, then barcode.

    Selection is independent of PAM50 labels, expression values, outcomes, and
    any downstream analysis. The caller must resolve or exclude ambiguous
    sample-type and patient-barcode evidence before calling this function.
    """
    grouped: dict[str, list[SampleAuditRecord]] = defaultdict(list)
    for record in samples:
        if record.clinical_match_status != "matched" or record.patient_id is None:
            continue
        if record.sample_type_priority is None or record.sample_type_evidence_status not in {
            "consistent",
            "partial_but_consistent",
        }:
            raise ValueError(
                f"Cannot rank sample {record.sample_id}: sample-type evidence is ambiguous"
            )
        if record.patient_barcode_validation_status == "mismatch":
            raise ValueError(
                f"Cannot group sample {record.sample_id}: clinical _PATIENT conflicts with barcode"
            )
        grouped[record.patient_id].append(record)

    selected: dict[str, SampleAuditRecord] = {}
    for patient_id, group in grouped.items():
        selected[patient_id] = min(
            group,
            key=lambda record: (
                record.sample_type_priority,
                record.sample_id.casefold(),
                record.sample_id,
            ),
        )
    return dict(sorted(selected.items()))


def build_classification_cohort(
    expression_sample_ids: Sequence[str],
    clinical_rows: Iterable[Mapping[str, str | None]],
) -> CohortBuildResult:
    """Build a primary PAM50 cohort and a disposition for every expression ID.

    One sample is selected per patient *before* checking PAM50 availability.
    This prevents a non-primary sample from replacing a selected primary tumor
    simply because its label is present. A missing or unmapped label on the
    selected sample is reported as an exclusion; there is no label-based
    fallback to another sample.
    """
    expression_ids = [sample_id.strip() for sample_id in expression_sample_ids]
    if any(not sample_id for sample_id in expression_ids):
        raise ValueError("Expression sample IDs must be nonblank")
    duplicates = [sample_id for sample_id, count in Counter(expression_ids).items() if count > 1]
    if duplicates:
        examples = ", ".join(sorted(duplicates)[:5])
        raise ValueError(f"Expression header contains duplicate sample IDs: {examples}")

    clinical_by_sample: dict[str, list[Mapping[str, str | None]]] = defaultdict(list)
    clinical_ids: set[str] = set()
    clinical_row_count = 0
    for row in clinical_rows:
        clinical_row_count += 1
        sample_id = (row.get("sampleID") or "").strip()
        if sample_id:
            clinical_ids.add(sample_id)
            clinical_by_sample[sample_id].append(row)

    expression_id_set = set(expression_ids)
    duplicate_clinical_ids = sorted(
        sample_id for sample_id, matched_rows in clinical_by_sample.items() if len(matched_rows) > 1
    )
    records: list[SampleAuditRecord] = []
    patient_counts: Counter[str] = Counter()

    for sample_id in expression_ids:
        patient_id = derive_patient_id(sample_id)
        if patient_id:
            patient_counts[patient_id] += 1

        matched_rows = clinical_by_sample.get(sample_id, [])
        if not matched_rows:
            record = _make_record(sample_id, patient_id, None, "no_match")
            record.sample_retention_reason = "no_matching_clinical_metadata"
        elif len(matched_rows) > 1:
            record = _make_record(sample_id, patient_id, None, "duplicate_clinical_sample_id")
            record.sample_retention_reason = "duplicate_clinical_sample_id_requires_review"
        else:
            record = _make_record(sample_id, patient_id, matched_rows[0], "matched")
            if patient_id is None:
                record.sample_retention_reason = "unparseable_tcga_patient_id"
        records.append(record)

    groups: dict[str, list[SampleAuditRecord]] = defaultdict(list)
    for record in records:
        if record.clinical_match_status == "matched" and record.patient_id is not None:
            groups[record.patient_id].append(record)

    selectable_records: list[SampleAuditRecord] = []
    for _patient_id, group in groups.items():
        if any(record.patient_barcode_validation_status == "mismatch" for record in group):
            for record in group:
                record.sample_retention_reason = "patient_identifier_conflict_requires_review"
            continue
        if any(
            record.sample_type_priority is None
            or record.sample_type_evidence_status not in {"consistent", "partial_but_consistent"}
            for record in group
        ):
            for record in group:
                record.sample_retention_reason = "sample_type_conflict_requires_review"
            continue
        selectable_records.extend(group)

    selected_by_patient = select_one_sample_per_patient(selectable_records)
    for patient_id, group in groups.items():
        if any(record.sample_retention_reason.endswith("requires_review") for record in group):
            continue
        selected = selected_by_patient[patient_id]
        for record in group:
            record.selected_sample_id_for_patient = selected.sample_id
            if record is selected:
                record.is_selected_patient_sample = True
                if record.pam50_label_status == "missing":
                    record.sample_retention_reason = "excluded_missing_pam50"
                elif record.pam50_label_status == "unmapped":
                    record.sample_retention_reason = "unmapped_pam50_label_requires_review"
                else:
                    record.is_in_classification_cohort = True
                    record.sample_retention_reason = f"retained_{record.sample_type_category}"
            elif record.sample_type_priority == selected.sample_type_priority:
                record.sample_retention_reason = "not_selected_lexicographic_tie_break"
            else:
                record.sample_retention_reason = "not_selected_lower_priority_sample_type"

    cohort = sorted(
        (record for record in records if record.is_in_classification_cohort),
        key=lambda record: (record.patient_id or "", record.sample_id),
    )
    patient_distribution = Counter(patient_counts.values())
    return CohortBuildResult(
        sample_records=records,
        cohort=cohort,
        expression_sample_count=len(expression_ids),
        clinical_row_count=clinical_row_count,
        clinical_unique_sample_count=len(clinical_ids),
        clinical_only_sample_count=len(clinical_ids - expression_id_set),
        duplicate_clinical_sample_ids=duplicate_clinical_ids,
        multiple_sample_patient_count=sum(
            patient_count
            for sample_count, patient_count in patient_distribution.items()
            if sample_count > 1
        ),
        patient_sample_count_distribution=dict(sorted(patient_distribution.items())),
    )
