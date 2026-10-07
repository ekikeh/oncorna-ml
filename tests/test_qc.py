from __future__ import annotations

import csv
import gzip

import numpy as np
import pytest

from oncorna.qc import (
    GENE_SUMMARY_FIELDS,
    SAMPLE_METRIC_FIELDS,
    ExpressionMatrix,
    build_gene_filter_summary,
    calculate_sample_metrics,
    read_classification_cohort,
    read_expression_matrix,
    run_expression_qc,
    tukey_fence_flags,
    write_tsv,
)


def test_sample_metrics_calculate_expression_statistics_and_deterministic_exclusions() -> None:
    values = np.asarray(
        [
            [0.0, 0.05, 2.0],
            [1.0, 2.0, np.nan],
            [4.0, np.inf, 3.0],
        ],
        dtype=np.float32,
    )

    rows, lower, upper = calculate_sample_metrics(
        values,
        ["S1", "S2", "S3"],
        selected_expression_threshold=1.0,
        near_zero_expression_max=0.1,
        maximum_missing_fraction=0.05,
        exclude_any_infinite_values=True,
        exclude_no_finite_values=True,
        tukey_iqr_multiplier=1.5,
    )

    first, second, third = rows
    assert first["aggregate_expression_signal"] == pytest.approx(5.0)
    assert first["median_expression"] == pytest.approx(1.0)
    assert first["mean_expression"] == pytest.approx(5.0 / 3.0)
    assert first["median_absolute_deviation"] == pytest.approx(1.0)
    assert first["n_genes_gt_0"] == 2
    assert first["n_genes_gt_1"] == 1
    assert first["n_genes_gt_2"] == 1
    assert first["n_zero_expression"] == 1
    assert first["n_near_zero_expression"] == 1
    assert first["sample_included_in_qc_cohort"] is True
    assert first["sample_exclusion_reason"] == "retained_for_qc"
    assert first["tukey_outlier_aggregate_signal"] is False
    assert lower == upper == pytest.approx(5.0)

    assert second["n_infinite_values"] == 1
    assert second["sample_included_in_qc_cohort"] is False
    assert second["sample_exclusion_reason"] == "excluded_infinite_expression_values"
    assert second["tukey_outlier_aggregate_signal"] is None
    assert third["missing_fraction"] == pytest.approx(1 / 3)
    assert third["sample_included_in_qc_cohort"] is False
    assert third["sample_exclusion_reason"] == "excluded_excess_missingness"


def test_tukey_fence_flags_are_strict_and_deterministic() -> None:
    values = np.asarray([10.0, 11.0, 12.0, 11.0, 100.0])

    first_flags, lower, upper = tukey_fence_flags(values, multiplier=1.5)
    second_flags, lower_again, upper_again = tukey_fence_flags(values, multiplier=1.5)

    assert first_flags.tolist() == [False, False, False, False, True]
    assert (lower, upper) == pytest.approx((9.5, 13.5))
    assert (lower_again, upper_again) == pytest.approx((lower, upper))
    assert np.array_equal(first_flags, second_flags)


def test_gene_prevalence_filter_is_deterministic_and_uses_strict_thresholds() -> None:
    gene_ids = ["g1", "g2", "g3"]
    values = np.asarray(
        [
            [0.0, 2.0, 2.0, 2.0],
            [2.0, 2.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0],
        ],
        dtype=np.float32,
    )
    included = [True, True, True, True]
    kwargs = {
        "sensitivity_thresholds": (0.0, 1.0, 2.0),
        "selected_expression_threshold": 1.0,
        "minimum_sample_fraction": 0.5,
    }

    first = build_gene_filter_summary(gene_ids, values, included, **kwargs)
    second = build_gene_filter_summary(gene_ids, values, included, **kwargs)

    assert first.keep_mask.tolist() == [True, True, False]
    assert first.sensitivity_counts == {0.0: 2, 1.0: 2, 2.0: 0}
    assert first.required_sample_count == 2
    assert np.array_equal(first.keep_mask, second.keep_mask)
    assert first.rows == second.rows
    assert first.rows[0]["n_samples_gt_2"] == 0


def test_gene_filter_has_no_pam50_input_and_class_changes_do_not_change_genes() -> None:
    sample_ids = ["S1", "S2", "S3", "S4"]
    values = np.asarray(
        [[0.0, 2.0, 2.0, 2.0], [2.0, 2.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0]],
        dtype=np.float32,
    )
    expression = ExpressionMatrix(
        gene_ids=("g1", "g2", "g3"),
        sample_ids=tuple(sample_ids),
        values=values,
        source_sample_count=4,
        noncohort_source_sample_count=0,
    )
    labels_a = ["Luminal A", "Luminal B", "Basal-like", "Normal-like"]
    labels_b = list(reversed(labels_a))

    def cohort_with_labels(labels: list[str]) -> list[dict[str, str]]:
        return [
            {
                "sample_id": sample_id,
                "patient_id": f"P{index}",
                "pam50_original_label": label,
                "pam50_normalized_label": label,
                "sample_type": "Primary Tumor",
            }
            for index, (sample_id, label) in enumerate(zip(sample_ids, labels, strict=True))
        ]

    result_a = run_expression_qc(
        expression,
        cohort_with_labels(labels_a),
        near_zero_expression_max=0.1,
        sensitivity_thresholds=(0.0, 1.0, 2.0),
        selected_expression_threshold=1.0,
        minimum_sample_fraction=0.5,
        maximum_missing_fraction=0.05,
        exclude_any_infinite_values=True,
        exclude_no_finite_values=True,
        tukey_iqr_multiplier=1.5,
    )
    result_b = run_expression_qc(
        expression,
        cohort_with_labels(labels_b),
        near_zero_expression_max=0.1,
        sensitivity_thresholds=(0.0, 1.0, 2.0),
        selected_expression_threshold=1.0,
        minimum_sample_fraction=0.5,
        maximum_missing_fraction=0.05,
        exclude_any_infinite_values=True,
        exclude_no_finite_values=True,
        tukey_iqr_multiplier=1.5,
    )

    assert result_a.gene_filter.keep_mask.tolist() == result_b.gene_filter.keep_mask.tolist()
    assert result_a.gene_filter.sensitivity_counts == result_b.gene_filter.sensitivity_counts


def test_cohort_loader_rejects_duplicate_patient_ids(tmp_path) -> None:
    cohort_path = tmp_path / "cohort.tsv"
    cohort_path.write_text(
        "sample_id\tpatient_id\tpam50_original_label\tpam50_normalized_label\tsample_type\n"
        "S1\tP1\tLumA\tLuminal A\tPrimary Tumor\n"
        "S2\tP1\tLumB\tLuminal B\tPrimary Tumor\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="duplicate patient IDs"):
        read_classification_cohort(cohort_path)


def test_expression_loader_reads_only_requested_cohort_columns(tmp_path) -> None:
    matrix_path = tmp_path / "tiny.tsv.gz"
    with gzip.open(matrix_path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample", "S1", "S2", "outside"])
        writer.writerow(["gene-a", "1.0", "2.0", "9999.0"])
        writer.writerow(["gene-b", "0.0", "3.0", "9999.0"])

    result = read_expression_matrix(matrix_path, ["S2", "S1"])

    assert result.sample_ids == ("S2", "S1")
    assert result.source_sample_count == 3
    assert result.noncohort_source_sample_count == 1
    assert result.gene_ids == ("gene-a", "gene-b")
    assert result.values.tolist() == [[2.0, 1.0], [3.0, 0.0]]


def test_qc_output_schema_is_stable(tmp_path) -> None:
    sample_path = tmp_path / "qc_sample_metrics.tsv"
    gene_path = tmp_path / "qc_gene_summary.tsv"
    write_tsv(sample_path, [{field: "x" for field in SAMPLE_METRIC_FIELDS}], SAMPLE_METRIC_FIELDS)
    write_tsv(gene_path, [{field: "x" for field in GENE_SUMMARY_FIELDS}], GENE_SUMMARY_FIELDS)

    with sample_path.open(newline="", encoding="utf-8") as handle:
        sample_header = next(csv.reader(handle, delimiter="\t"))
    with gene_path.open(newline="", encoding="utf-8") as handle:
        gene_header = next(csv.reader(handle, delimiter="\t"))

    assert sample_header == list(SAMPLE_METRIC_FIELDS)
    assert gene_header == list(GENE_SUMMARY_FIELDS)
