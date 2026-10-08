from __future__ import annotations

import csv
import gzip
import json
import math
from pathlib import Path

import pytest

from oncorna.ml_matrix import (
    GENE_SUMMARY_REQUIRED_FIELDS,
    _load_gene_filter_results,
    build_ml_matrix,
)

COHORT_FIELDS = (
    "sample_id",
    "patient_id",
    "pam50_original_label",
    "pam50_normalized_label",
    "sample_type",
)
SAMPLE_METRIC_FIELDS = (
    "sample_id",
    "sample_included_in_qc_cohort",
    "sample_exclusion_reason",
    "tukey_outlier_aggregate_signal",
    "n_missing_values",
    "n_infinite_values",
)


def _write_tsv(path: Path, fields: tuple[str, ...], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def _write_expression(path: Path, header: list[str], rows: list[list[str]]) -> None:
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def _gene_summary_row(
    gene_id: str,
    n_above: int,
    *,
    sample_count: int = 2,
    required: int = 1,
    threshold: float = 1.0,
) -> dict[str, str]:
    return {
        "gene_id": gene_id,
        "n_samples_in_qc_cohort": str(sample_count),
        "n_samples_above_selected_threshold": str(n_above),
        "fraction_samples_above_selected_threshold": str(n_above / sample_count),
        "selected_filter_expression_gt": str(threshold),
        "required_samples_for_selected_filter": str(required),
        "retained_by_prevalence_filter": str(n_above >= required).lower(),
    }


@pytest.fixture
def mini_project(tmp_path: Path) -> dict[str, Path]:
    cohort_path = tmp_path / "data/processed/classification_cohort.tsv"
    _write_tsv(
        cohort_path,
        COHORT_FIELDS,
        [
            {
                "sample_id": "S1",
                "patient_id": "P1",
                "pam50_original_label": "Normal",
                "pam50_normalized_label": "Normal-like",
                "sample_type": "Primary Tumor",
            },
            {
                "sample_id": "S2",
                "patient_id": "P2",
                "pam50_original_label": "LumA",
                "pam50_normalized_label": "Luminal A",
                "sample_type": "Primary Tumor",
            },
        ],
    )

    sample_metrics_path = tmp_path / "data/processed/qc_sample_metrics.tsv"
    _write_tsv(
        sample_metrics_path,
        SAMPLE_METRIC_FIELDS,
        [
            {
                "sample_id": "S1",
                "sample_included_in_qc_cohort": "true",
                "sample_exclusion_reason": "retained_for_qc",
                "tukey_outlier_aggregate_signal": "false",
                "n_missing_values": "0",
                "n_infinite_values": "0",
            },
            {
                "sample_id": "S2",
                "sample_included_in_qc_cohort": "true",
                "sample_exclusion_reason": "retained_for_qc",
                "tukey_outlier_aggregate_signal": "true",
                "n_missing_values": "0",
                "n_infinite_values": "0",
            },
        ],
    )

    gene_summary_path = tmp_path / "data/processed/qc_gene_summary.tsv"
    _write_tsv(
        gene_summary_path,
        GENE_SUMMARY_REQUIRED_FIELDS,
        [
            _gene_summary_row("G_SECOND", 1),
            _gene_summary_row("G_FIRST", 1),
            _gene_summary_row("G_DROP", 0),
        ],
    )

    expression_path = tmp_path / "data/raw/TCGA-BRCA_HiSeqV2.tsv.gz"
    expression_path.parent.mkdir(parents=True, exist_ok=True)
    # Source order differs from cohort order; OUTSIDE is not an approved sample.
    _write_expression(
        expression_path,
        ["sample", "OUTSIDE", "S2", "S1"],
        [
            ["G_SECOND", "999.0", "2.50", "0.25"],
            ["G_FIRST", "999.0", "3.40", "1.0"],
            ["G_DROP", "999.0", "1.0", "0.0"],
        ],
    )

    config_path = tmp_path / "configs/default.yaml"
    config_path.parent.mkdir(parents=True, exist_ok=True)
    config_path.write_text("test-config: phase-1e\n", encoding="utf-8")

    return {
        "root": tmp_path,
        "cohort": cohort_path,
        "sample_metrics": sample_metrics_path,
        "gene_summary": gene_summary_path,
        "expression": expression_path,
        "config": config_path,
        "matrix": tmp_path / "data/processed/ml_expression_matrix.tsv.gz",
        "gene_list": tmp_path / "data/processed/ml_gene_list.tsv",
        "manifest": tmp_path / "data/processed/ml_matrix_manifest.json",
    }


def _build(paths: dict[str, Path], **overrides: object):
    arguments: dict[str, object] = {
        "expression_path": paths["expression"],
        "cohort_path": paths["cohort"],
        "gene_summary_path": paths["gene_summary"],
        "sample_metrics_path": paths["sample_metrics"],
        "config_path": paths["config"],
        "output_path": paths["matrix"],
        "gene_list_path": paths["gene_list"],
        "manifest_path": paths["manifest"],
        "project_root": paths["root"],
        "expression_scale": "log2(normalized_count + 1)",
        "expression_threshold": 1.0,
        "minimum_sample_fraction": 0.5,
        "expected_sample_count": 2,
        "expected_gene_count": 2,
        "expected_source_gene_count": 3,
        "expression_transform": "none; copy selected source values verbatim",
        "config_snapshot": {"qc": {"gene_filter": {"expression_gt": 1.0}}},
    }
    arguments.update(overrides)
    return build_ml_matrix(**arguments)  # type: ignore[arg-type]


def test_matrix_uses_only_cohort_samples_and_preserves_values_and_order(
    mini_project: dict[str, Path],
) -> None:
    result = _build(mini_project)

    assert result.sample_ids == ("S1", "S2")
    assert result.gene_ids == ("G_SECOND", "G_FIRST")
    assert result.source_gene_count == 3
    assert result.source_sample_count == 3
    assert result.noncohort_source_sample_count == 1
    assert result.required_sample_count == 1
    assert result.flagged_outlier_sample_ids == ("S2",)

    with gzip.open(result.matrix_path, "rt", encoding="utf-8", newline="") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    assert rows == [
        ["gene_id", "S1", "S2"],
        ["G_SECOND", "0.25", "2.50"],
        ["G_FIRST", "1.0", "3.40"],
    ]
    assert all("OUTSIDE" not in field for row in rows for field in row)

    # These are source tokens, not a second log transform of the processed values.
    assert rows[1][1:] == ["0.25", "2.50"]
    assert float(rows[1][2]) != pytest.approx(math.log2(2.50 + 1.0))

    with mini_project["gene_list"].open("r", encoding="utf-8", newline="") as handle:
        gene_rows = list(csv.DictReader(handle, delimiter="\t"))
    assert [row["gene_id"] for row in gene_rows] == ["G_SECOND", "G_FIRST"]
    assert "G_DROP" not in {row["gene_id"] for row in gene_rows}


def test_manifest_records_dimensions_inputs_configuration_and_orders(
    mini_project: dict[str, Path],
) -> None:
    _build(mini_project)
    manifest = json.loads(mini_project["manifest"].read_text(encoding="utf-8"))

    assert manifest["expression_units"] == "log2(normalized_count + 1)"
    assert manifest["source_matrix"]["gene_rows"] == 3
    assert manifest["source_matrix"]["sample_columns"] == 3
    assert manifest["source_matrix"]["noncohort_sample_columns_omitted"] == 1
    assert manifest["matrix"]["sample_order_source"] == "classification_cohort.tsv row order"
    assert manifest["matrix"]["additional_transformation"].startswith("none;")
    assert manifest["matrix"]["shape"] == {
        "gene_rows": 2,
        "sample_columns": 2,
        "total_columns_including_gene_id": 3,
    }
    assert manifest["sample_order"] == ["S1", "S2"]
    assert manifest["gene_order"] == ["G_SECOND", "G_FIRST"]
    assert manifest["tukey_outlier_samples_flagged_but_retained"] == ["S2"]
    assert manifest["clinical_or_pam50_columns_in_expression_matrix"] is False
    assert manifest["splits_or_modeling_performed"] is False
    assert {item["path"] for item in manifest["inputs"]} == {
        "data/raw/TCGA-BRCA_HiSeqV2.tsv.gz",
        "data/processed/classification_cohort.tsv",
        "data/processed/qc_gene_summary.tsv",
        "data/processed/qc_sample_metrics.tsv",
        "configs/default.yaml",
    }
    assert manifest["configuration"]["qc"]["gene_filter"]["expression_gt"] == 1.0


def test_configured_phase_1d_prevalence_boundary_is_169_of_844(tmp_path: Path) -> None:
    path = tmp_path / "qc_gene_summary.tsv"
    rows = [
        _gene_summary_row("edge_kept", 169, sample_count=844, required=169),
        _gene_summary_row("below_edge", 168, sample_count=844, required=169),
    ]
    _write_tsv(path, GENE_SUMMARY_REQUIRED_FIELDS, rows)

    all_rows, selected_rows, required = _load_gene_filter_results(
        path,
        expected_sample_count=844,
        expected_source_gene_count=2,
        expression_threshold=1.0,
        minimum_sample_fraction=0.20,
        expected_gene_count=1,
    )

    assert required == 169
    assert len(all_rows) == 2
    assert [row["gene_id"] for row in selected_rows] == ["edge_kept"]


def test_duplicate_cohort_sample_ids_are_rejected(mini_project: dict[str, Path]) -> None:
    _write_tsv(
        mini_project["cohort"],
        COHORT_FIELDS,
        [
            {
                "sample_id": "S1",
                "patient_id": "P1",
                "pam50_original_label": "LumA",
                "pam50_normalized_label": "Luminal A",
                "sample_type": "Primary Tumor",
            },
            {
                "sample_id": "S1",
                "patient_id": "P2",
                "pam50_original_label": "Normal",
                "pam50_normalized_label": "Normal-like",
                "sample_type": "Primary Tumor",
            },
        ],
    )

    with pytest.raises(ValueError, match="duplicate sample IDs"):
        _build(mini_project)


def test_duplicate_source_sample_ids_are_rejected(mini_project: dict[str, Path]) -> None:
    _write_expression(
        mini_project["expression"],
        ["sample", "OUTSIDE", "S1", "S1"],
        [["G_SECOND", "999.0", "0.25", "2.50"]],
    )

    with pytest.raises(ValueError, match="duplicate column identifiers"):
        _build(mini_project)


def test_duplicate_source_gene_ids_are_rejected(mini_project: dict[str, Path]) -> None:
    _write_expression(
        mini_project["expression"],
        ["sample", "OUTSIDE", "S2", "S1"],
        [
            ["G_SECOND", "999.0", "2.50", "0.25"],
            ["G_SECOND", "999.0", "2.50", "0.25"],
        ],
    )

    with pytest.raises(ValueError, match="Duplicate source gene ID"):
        _build(mini_project)


def test_duplicate_gene_ids_in_phase_1d_summary_are_rejected(
    mini_project: dict[str, Path],
) -> None:
    _write_tsv(
        mini_project["gene_summary"],
        GENE_SUMMARY_REQUIRED_FIELDS,
        [
            _gene_summary_row("G_DUP", 1),
            _gene_summary_row("G_DUP", 1),
            _gene_summary_row("G_DROP", 0),
        ],
    )

    with pytest.raises(ValueError, match="Duplicate gene_id"):
        _build(mini_project)


def test_non_none_transformation_is_rejected(mini_project: dict[str, Path]) -> None:
    with pytest.raises(ValueError, match="does not apply a second expression transformation"):
        _build(mini_project, expression_transform="log2 again")


def test_raw_count_units_are_rejected(mini_project: dict[str, Path]) -> None:
    with pytest.raises(
        ValueError, match=r"requires expression units log2\(normalized_count \+ 1\)"
    ):
        _build(mini_project, expression_scale="raw counts")


def test_outputs_must_remain_under_git_ignored_processed_directory(
    mini_project: dict[str, Path],
) -> None:
    with pytest.raises(ValueError, match="Git-ignored data/processed"):
        _build(mini_project, output_path=mini_project["root"] / "ml_expression_matrix.tsv.gz")
