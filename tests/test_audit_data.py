from __future__ import annotations

import csv
import gzip
from pathlib import Path

from scripts import audit_data


def _write_tsv(path: Path, rows: list[list[str]], *, compressed: bool = False) -> None:
    opener = gzip.open if compressed else Path.open
    with opener(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerows(rows)


def test_parse_tcga_patient_id_uses_first_three_barcode_fields() -> None:
    assert audit_data.parse_tcga_patient_id("TCGA-AR-A5QQ-01") == "TCGA-AR-A5QQ"
    assert audit_data.parse_tcga_patient_id(" TCGA-AR-A5QQ-11A ") == "TCGA-AR-A5QQ"
    assert audit_data.parse_tcga_patient_id("SYN-AR-A5QQ-01") is None
    assert audit_data.parse_tcga_patient_id("TCGA-A-A5QQ-01") is None


def test_classify_gene_identifier_describes_visible_patterns_only() -> None:
    assert audit_data.classify_gene_identifier("BRCA1") == "Gene-symbol-like identifier"
    assert audit_data.classify_gene_identifier("ENSG000001234.5") == "Ensembl gene ID"
    assert audit_data.classify_gene_identifier("?|12345") == "Xena fallback with numeric ID"
    assert audit_data.classify_gene_identifier("") == "Blank or missing identifier"


def test_key_column_detection_avoids_substring_false_positives() -> None:
    columns = [
        "person_neoplasm_cancer_status",
        "ER_Status_nature2012",
        "er_level_cell_percentage_category",
        "metastatic_breast_carcinoma_estrogen_receptor_detection_mthd_txt",
        "positive_finding_estrogen_receptor_other_measurement_scale_text",
        "Age_at_Initial_Pathologic_Diagnosis_nature2012",
        "AJCC_Stage_nature2012",
        "days_to_birth",
        "PAM50Call_RNAseq",
    ]

    categories = audit_data._key_columns(columns)

    assert "person_neoplasm_cancer_status" not in categories["ER"]
    assert "ER_Status_nature2012" in categories["ER"]
    assert "er_level_cell_percentage_category" in categories["ER"]
    assert "metastatic_breast_carcinoma_estrogen_receptor_detection_mthd_txt" in categories["ER"]
    assert "positive_finding_estrogen_receptor_other_measurement_scale_text" in categories["ER"]
    assert "er_level_cell_percentage_category" not in categories["Age"]
    assert "AJCC_Stage_nature2012" not in categories["Age"]
    assert "Age_at_Initial_Pathologic_Diagnosis_nature2012" in categories["Age"]
    assert "days_to_birth" in categories["Age"]
    assert "PAM50Call_RNAseq" in categories["PAM50 subtype"]


def test_expression_inspection_counts_dimensions_and_duplicate_genes(tmp_path: Path) -> None:
    expression_path = tmp_path / "expression.tsv.gz"
    _write_tsv(
        expression_path,
        [
            ["sample", "TCGA-AA-1234-01", "TCGA-AA-1234-11"],
            ["BRCA1", "1.0", "2.0"],
            ["ENSG000001234.5", "3.0", "4.0"],
            ["BRCA1", "5.0", "6.0"],
        ],
        compressed=True,
    )

    result = audit_data.inspect_expression(expression_path)

    assert result.identifier_column == "sample"
    assert result.sample_columns == 2
    assert result.gene_rows == 3
    assert result.nonempty_gene_ids == 3
    assert result.unique_gene_ids == 2
    assert result.duplicate_gene_ids == {"BRCA1": 2}
    assert result.duplicate_gene_rows == 1
    assert result.identifier_types == {
        "Ensembl gene ID": 1,
        "Gene-symbol-like identifier": 2,
    }
    assert result.malformed_rows == 0


def test_run_audit_reports_matching_missingness_labels_and_patient_multiplicity(
    tmp_path: Path,
) -> None:
    expression_path = tmp_path / "expression.tsv.gz"
    clinical_path = tmp_path / "clinical.tsv"
    report_path = tmp_path / "reports" / "data_audit.md"
    expression_ids = [
        "TCGA-AA-1234-01",
        "TCGA-AA-1234-11",
        "TCGA-BB-5678-01",
    ]
    _write_tsv(
        expression_path,
        [
            ["sample", *expression_ids],
            ["BRCA1", "1.0", "2.0", "3.0"],
            ["ENSG000001234.5", "4.0", "5.0", "6.0"],
        ],
        compressed=True,
    )
    _write_tsv(
        clinical_path,
        [
            [
                "sampleID",
                "patient_id",
                "bcr_patient_barcode",
                "_PATIENT",
                "PAM50Call_RNAseq",
                "PAM50_mRNA_nature2012",
                "ER_Status_nature2012",
                "person_neoplasm_cancer_status",
                "Age_at_Initial_Pathologic_Diagnosis_nature2012",
                "AJCC_Stage_nature2012",
            ],
            [
                "TCGA-AA-1234-01",
                "1234",
                "TCGA-AA-1234",
                "TCGA-AA-1234",
                "LumA",
                "Luminal A",
                "Positive",
                "Alive",
                "55",
                "Stage II",
            ],
            [
                "TCGA-BB-5678-01",
                "5678",
                "TCGA-BB-5678",
                "TCGA-BB-5678",
                "Basal",
                "Basal-like",
                "NA",
                "Dead",
                "68",
                "Stage I",
            ],
            [
                "TCGA-CC-9999-01",
                "9999",
                "TCGA-CC-9999",
                "TCGA-CC-9999",
                " lumA ",
                "",
                "Negative",
                "Alive",
                "",
                "Stage III",
            ],
            [
                "TCGA-DD-0001-01",
                "0001",
                "TCGA-DD-0001",
                "TCGA-DD-0001",
                "",
                "",
                "",
                "Alive",
                "60",
                "Stage I",
            ],
        ],
    )

    result = audit_data.run_audit(
        expression_path=expression_path,
        clinical_path=clinical_path,
        report_path=report_path,
        audit_date="2026-10-07",
    )

    assert result.matching.intersection_ids == ["TCGA-AA-1234-01", "TCGA-BB-5678-01"]
    assert result.matching.expression_only_ids == ["TCGA-AA-1234-11"]
    assert result.matching.clinical_only_ids == ["TCGA-CC-9999-01", "TCGA-DD-0001-01"]
    assert result.patients.patients_with_multiple_samples == 1
    assert result.patients.samples_per_patient_distribution == {1: 1, 2: 1}
    assert result.clinical.patient_suffix_comparisons == 4
    assert result.clinical.patient_suffix_matches == 4
    assert result.clinical.pam50_counts == {" lumA ": 1, "Basal": 1, "LumA": 1}
    assert result.clinical.pam50_missing == 1
    assert result.clinical.pam50_case_or_spacing_variants == [[" lumA ", "LumA"]]
    assert result.clinical.missing_counts["ER_Status_nature2012"] == 2
    assert "person_neoplasm_cancer_status" not in result.clinical.key_columns["ER"]
    assert "AJCC_Stage_nature2012" not in result.clinical.key_columns["Age"]
    assert report_path.is_file()
    assert "66.7%" in result.markdown
    assert "50.0%" in result.markdown
    assert "`LumA`: 1" in result.markdown
    assert "` lumA `: 1" in result.markdown
    assert "Ensembl-style identifiers detected by this pattern check: **1**" in result.markdown
    assert "`patient_id` `1234` matches the final segment" in result.markdown
    assert "Patients with more than one expression sample: **1**" in result.markdown
