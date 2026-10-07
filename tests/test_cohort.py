from __future__ import annotations

from collections import Counter

import pytest

from oncorna.cohort import (
    build_classification_cohort,
    classify_sample_type,
    derive_patient_id,
    normalize_pam50,
    parse_sample_type_code,
)


def _clinical_row(
    sample_id: str,
    label: str,
    sample_type: str,
    sample_type_id: str,
) -> dict[str, str]:
    patient_id = derive_patient_id(sample_id)
    assert patient_id is not None
    sample_type_code = parse_sample_type_code(sample_id)
    assert sample_type_code is not None
    return {
        "sampleID": sample_id,
        "PAM50Call_RNAseq": label,
        "sample_type": sample_type,
        "sample_type_id": sample_type_id,
        "_PATIENT": patient_id,
        "patient_id": patient_id.rsplit("-", maxsplit=1)[-1],
        "bcr_patient_barcode": patient_id,
        "bcr_sample_barcode": f"{patient_id}-{sample_type_code}A",
    }


def test_patient_and_sample_type_codes_parse_from_tcga_barcode() -> None:
    barcode = "TCGA-AR-A5QQ-01A-11D-A22X-01"

    assert derive_patient_id(barcode) == "TCGA-AR-A5QQ"
    assert parse_sample_type_code(barcode) == "01"
    assert derive_patient_id("TCGA-A-A5QQ-01") is None
    assert parse_sample_type_code("TCGA-AR-A5QQ") is None
    assert derive_patient_id("SYN-AR-A5QQ-01") is None


def test_normalize_pam50_maps_all_five_labels_without_mutating_raw_values() -> None:
    expected = {
        "LumA": "Luminal A",
        "LumB": "Luminal B",
        "Basal": "Basal-like",
        "Her2": "HER2-enriched",
        "Normal": "Normal-like",
    }

    for original, normalized in expected.items():
        assert normalize_pam50(original) == normalized
    assert normalize_pam50("  LumA  ") == "Luminal A"
    assert normalize_pam50("") is None
    assert normalize_pam50("not-a-source-label") is None


def test_sample_type_classification_checks_barcode_and_clinical_evidence() -> None:
    primary = classify_sample_type("TCGA-AA-1234-01", "Primary Tumor", "01")
    metastatic = classify_sample_type("TCGA-BB-5678-06", "Metastatic", "06")
    normal = classify_sample_type("TCGA-CC-9012-11", "Solid Tissue Normal", "11")
    conflict = classify_sample_type("TCGA-AA-1234-01", "Solid Tissue Normal", "11")

    assert (primary.category, primary.priority, primary.evidence_status) == (
        "primary_tumor",
        0,
        "consistent",
    )
    assert metastatic.category == "metastatic"
    assert normal.category == "solid_tissue_normal"
    assert conflict.category == "ambiguous"
    assert conflict.priority is None
    assert conflict.evidence_status == "sample_type_sources_disagree"


def test_one_patient_is_retained_once_and_duplicate_resolution_is_deterministic() -> None:
    expression_ids = [
        "TCGA-AA-1234-11",
        "TCGA-AA-1234-01B",
        "TCGA-AA-1234-01A",
        "TCGA-BB-5678-01",
    ]
    clinical_rows = [
        _clinical_row("TCGA-AA-1234-11", "Normal", "Solid Tissue Normal", "11"),
        _clinical_row("TCGA-AA-1234-01B", "LumA", "Primary Tumor", "01"),
        _clinical_row("TCGA-AA-1234-01A", "Her2", "Primary Tumor", "01"),
        _clinical_row("TCGA-BB-5678-01", "Basal", "Primary Tumor", "01"),
    ]

    result = build_classification_cohort(expression_ids, clinical_rows)

    assert [row.sample_id for row in result.cohort] == [
        "TCGA-AA-1234-01A",
        "TCGA-BB-5678-01",
    ]
    assert len({row.patient_id for row in result.cohort}) == len(result.cohort) == 2
    assert result.cohort[0].pam50_original_label == "Her2"
    assert result.cohort[0].pam50_normalized_label == "HER2-enriched"
    assert result.cohort[0].sample_retention_reason == "retained_primary_tumor"

    reversed_result = build_classification_cohort(
        list(reversed(expression_ids)), list(reversed(clinical_rows))
    )
    assert {row.patient_id: row.sample_id for row in result.cohort} == {
        row.patient_id: row.sample_id for row in reversed_result.cohort
    }


def test_missing_pam50_on_preferred_primary_is_excluded_without_secondary_fallback() -> None:
    expression_ids = ["TCGA-CC-9012-01", "TCGA-CC-9012-11"]
    clinical_rows = [
        _clinical_row("TCGA-CC-9012-01", "NA", "Primary Tumor", "01"),
        _clinical_row("TCGA-CC-9012-11", "Normal", "Solid Tissue Normal", "11"),
    ]

    result = build_classification_cohort(expression_ids, clinical_rows)
    by_sample = {row.sample_id: row for row in result.sample_records}

    assert result.cohort == []
    assert result.selected_patient_sample_count == 1
    assert result.missing_pam50_matched_sample_count == 1
    assert result.missing_pam50_selected_sample_count == 1
    assert by_sample["TCGA-CC-9012-01"].sample_retention_reason == "excluded_missing_pam50"
    assert (
        by_sample["TCGA-CC-9012-11"].sample_retention_reason
        == "not_selected_lower_priority_sample_type"
    )


def test_derived_sample_patient_id_wins_over_truncated_legacy_patient_fields() -> None:
    sample_id = "TCGA-E9-A1NA-01"
    row = _clinical_row(sample_id, "LumB", "Primary Tumor", "01")
    row["patient_id"] = "A1"
    row["bcr_patient_barcode"] = "TCGA-E9-A1"
    row["bcr_sample_barcode"] = "TCGA-E9-A1-01A"

    result = build_classification_cohort([sample_id], [row])

    assert len(result.cohort) == 1
    record = result.cohort[0]
    assert record.patient_id == "TCGA-E9-A1NA"
    assert record.patient_barcode_validation_status == "match"
    assert record.clinical_patient_id_status == "mismatch"
    assert record.clinical_bcr_patient_barcode_status == "mismatch"
    assert record.bcr_sample_barcode_status == "mismatch"


def test_unmapped_label_and_sample_type_conflicts_are_explicit_not_silent() -> None:
    expression_ids = ["TCGA-DD-3456-01", "TCGA-EE-7890-01"]
    clinical_rows = [
        _clinical_row("TCGA-DD-3456-01", "Luminal A", "Primary Tumor", "01"),
        _clinical_row("TCGA-EE-7890-01", "LumA", "Solid Tissue Normal", "11"),
    ]

    result = build_classification_cohort(expression_ids, clinical_rows)
    by_sample = {row.sample_id: row for row in result.sample_records}

    assert result.cohort == []
    assert (
        by_sample["TCGA-DD-3456-01"].sample_retention_reason
        == "unmapped_pam50_label_requires_review"
    )
    assert (
        by_sample["TCGA-EE-7890-01"].sample_retention_reason
        == "sample_type_conflict_requires_review"
    )


def test_duplicate_expression_ids_fail_loudly() -> None:
    row = _clinical_row("TCGA-FF-1111-01", "LumA", "Primary Tumor", "01")

    with pytest.raises(ValueError, match="duplicate sample IDs"):
        build_classification_cohort(
            ["TCGA-FF-1111-01", "TCGA-FF-1111-01"],
            [row],
        )


def test_disposition_counts_are_sample_level_and_include_every_expression_id() -> None:
    expression_ids = ["TCGA-GG-2222-01", "TCGA-GG-2222-11", "TCGA-HH-3333-01"]
    clinical_rows = [
        _clinical_row("TCGA-GG-2222-01", "LumA", "Primary Tumor", "01"),
        _clinical_row("TCGA-GG-2222-11", "Normal", "Solid Tissue Normal", "11"),
        _clinical_row("TCGA-HH-3333-01", "", "Primary Tumor", "01"),
        _clinical_row("TCGA-II-4444-01", "Basal", "Primary Tumor", "01"),
    ]

    result = build_classification_cohort(expression_ids, clinical_rows)

    assert sum(result.disposition_counts.values()) == result.expression_sample_count == 3
    assert result.clinical_only_sample_count == 1
    assert result.matched_expression_count == 3
    assert result.multiple_sample_patient_count == 1
    assert Counter(row.sample_retention_reason for row in result.cohort) == {
        "retained_primary_tumor": 1
    }
