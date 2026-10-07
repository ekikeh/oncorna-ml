"""Expression and sample QC for the TCGA-BRCA HiSeqV2 matrix.

The source values are log2(normalized_count + 1), not raw integer counts.
This module computes expression-level descriptive statistics and prevalence
filters only; it does not implement count-based differential-expression QC.
"""

from __future__ import annotations

import csv
import gzip
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np

MISSING_TOKENS = {"", "na", "nan", "null", "none", "n/a", "--"}
SENSITIVITY_THRESHOLDS = (0.0, 1.0, 2.0)

SAMPLE_METRIC_FIELDS = (
    "sample_id",
    "patient_id",
    "pam50_original_label",
    "pam50_normalized_label",
    "sample_type",
    "n_genes_total",
    "n_finite_values",
    "n_missing_values",
    "n_infinite_values",
    "missing_fraction",
    "infinite_fraction",
    "nonfinite_fraction",
    "minimum_expression",
    "maximum_expression",
    "mean_expression",
    "median_expression",
    "q1_expression",
    "q3_expression",
    "iqr_expression",
    "median_absolute_deviation",
    "aggregate_expression_signal",
    "n_genes_gt_0",
    "n_genes_gt_1",
    "n_genes_gt_2",
    "n_genes_above_selected_threshold",
    "n_zero_expression",
    "zero_expression_fraction",
    "n_near_zero_expression",
    "near_zero_expression_fraction",
    "tukey_outlier_aggregate_signal",
    "aggregate_signal_lower_fence",
    "aggregate_signal_upper_fence",
    "sample_included_in_qc_cohort",
    "sample_exclusion_reason",
)

GENE_SUMMARY_FIELDS = (
    "gene_id",
    "n_samples_in_qc_cohort",
    "n_finite_values",
    "n_missing_values",
    "n_infinite_values",
    "minimum_expression",
    "maximum_expression",
    "mean_expression",
    "median_expression",
    "n_samples_gt_0",
    "fraction_samples_gt_0",
    "n_samples_gt_1",
    "fraction_samples_gt_1",
    "n_samples_gt_2",
    "fraction_samples_gt_2",
    "selected_filter_expression_gt",
    "n_samples_above_selected_threshold",
    "fraction_samples_above_selected_threshold",
    "required_samples_for_selected_filter",
    "retained_by_prevalence_filter",
)

COHORT_REQUIRED_FIELDS = (
    "sample_id",
    "patient_id",
    "pam50_original_label",
    "pam50_normalized_label",
    "sample_type",
)


@dataclass(frozen=True)
class ExpressionMatrix:
    """Selected expression matrix and source-header audit metadata."""

    gene_ids: tuple[str, ...]
    sample_ids: tuple[str, ...]
    values: np.ndarray
    source_sample_count: int
    noncohort_source_sample_count: int


@dataclass(frozen=True)
class GeneFilterResult:
    """Gene-level summary and a label-independent prevalence-filter mask."""

    rows: list[dict[str, Any]]
    sensitivity_counts: dict[float, int]
    required_sample_count: int
    keep_mask: np.ndarray


@dataclass(frozen=True)
class QCRunResult:
    """In-memory QC results; no transformed matrix is written by this module."""

    expression: ExpressionMatrix
    sample_metrics: list[dict[str, Any]]
    gene_filter: GeneFilterResult
    global_stats: dict[str, Any]
    outlier_lower_fence: float
    outlier_upper_fence: float


def read_classification_cohort(
    path: Path,
    *,
    expected_sample_count: int | None = None,
) -> list[dict[str, str]]:
    """Load and validate the frozen one-sample-per-patient cohort table."""
    path = Path(path)
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = set(reader.fieldnames or ())
        missing_fields = sorted(set(COHORT_REQUIRED_FIELDS) - fields)
        if missing_fields:
            raise ValueError(f"Cohort table is missing required columns: {missing_fields}")
        records = [
            {key: (value or "").strip() for key, value in row.items() if key is not None}
            for row in reader
        ]

    if not records:
        raise ValueError("Classification cohort contains no samples")
    sample_ids = [row["sample_id"] for row in records]
    patient_ids = [row["patient_id"] for row in records]
    if any(not value for value in sample_ids):
        raise ValueError("Classification cohort contains a blank sample_id")
    if any(not value for value in patient_ids):
        raise ValueError("Classification cohort contains a blank patient_id")
    if len(sample_ids) != len(set(sample_ids)):
        raise ValueError("Classification cohort contains duplicate sample IDs")
    if len(patient_ids) != len(set(patient_ids)):
        raise ValueError("Classification cohort contains duplicate patient IDs")
    if expected_sample_count is not None and len(records) != expected_sample_count:
        raise ValueError(
            "Classification cohort size does not match the frozen QC expectation: "
            f"expected {expected_sample_count}, found {len(records)}"
        )
    return records


def _parse_expression_value(raw: str, *, gene_id: str, sample_id: str) -> float:
    value = raw.strip()
    if value.casefold() in MISSING_TOKENS:
        return float("nan")
    try:
        return float(value)
    except ValueError as exc:
        raise ValueError(
            f"Non-numeric expression value for gene {gene_id!r}, sample {sample_id!r}: {raw!r}"
        ) from exc


def read_expression_matrix(
    path: Path,
    cohort_sample_ids: Sequence[str],
) -> ExpressionMatrix:
    """Read only cohort columns from a gzipped, gene-by-sample TSV matrix.

    The matrix is scanned twice to validate the full header/row structure before
    allocating a compact float32 array for the selected cohort columns. Source
    samples outside the supplied cohort are counted but never loaded.
    """
    path = Path(path)
    cohort_ids = [sample_id.strip() for sample_id in cohort_sample_ids]
    if not cohort_ids or any(not sample_id for sample_id in cohort_ids):
        raise ValueError("Cohort sample IDs must be a nonempty list of nonblank values")
    if len(cohort_ids) != len(set(cohort_ids)):
        raise ValueError("Cohort sample IDs contain duplicates")

    with gzip.open(path, "rt", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, None)
        if header is None or len(header) < 2:
            raise ValueError("Expression matrix must have a gene column and sample columns")
        header = [value.strip() for value in header]
        if len(header) != len(set(header)):
            raise ValueError("Expression matrix contains duplicate column IDs")
        sample_header = header[1:]
        index_by_sample = {sample_id: index + 1 for index, sample_id in enumerate(sample_header)}
        missing_ids = sorted(set(cohort_ids) - set(index_by_sample))
        if missing_ids:
            raise ValueError(
                f"Cohort sample IDs absent from expression matrix ({len(missing_ids)}): "
                f"{', '.join(missing_ids[:10])}"
            )
        selected_indices = [index_by_sample[sample_id] for sample_id in cohort_ids]
        gene_ids: list[str] = []
        seen_genes: set[str] = set()
        expected_width = len(header)
        for row_number, row in enumerate(reader, start=2):
            if len(row) != expected_width:
                raise ValueError(
                    f"Expression row {row_number} has {len(row)} fields; expected {expected_width}"
                )
            gene_id = row[0].strip()
            if not gene_id:
                raise ValueError(f"Expression row {row_number} has a blank gene identifier")
            if gene_id in seen_genes:
                raise ValueError(f"Expression matrix contains duplicate gene ID: {gene_id}")
            seen_genes.add(gene_id)
            gene_ids.append(gene_id)

    if not gene_ids:
        raise ValueError("Expression matrix contains no gene rows")

    matrix = np.empty((len(gene_ids), len(cohort_ids)), dtype=np.float32)
    with gzip.open(path, "rt", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        second_header = next(reader, None)
        if second_header is None or [value.strip() for value in second_header] != header:
            raise ValueError("Expression matrix header changed between validation and value passes")
        for row_index, row in enumerate(reader):
            if row_index >= len(gene_ids):
                raise ValueError(
                    "Expression matrix gained rows between validation and value passes"
                )
            gene_id = gene_ids[row_index]
            for column_index, source_index in enumerate(selected_indices):
                matrix[row_index, column_index] = _parse_expression_value(
                    row[source_index],
                    gene_id=gene_id,
                    sample_id=cohort_ids[column_index],
                )
        if row_index + 1 != len(gene_ids):
            raise ValueError("Expression matrix lost rows between validation and value passes")

    return ExpressionMatrix(
        gene_ids=tuple(gene_ids),
        sample_ids=tuple(cohort_ids),
        values=matrix,
        source_sample_count=len(sample_header),
        noncohort_source_sample_count=len(sample_header) - len(cohort_ids),
    )


def calculate_global_expression_stats(
    values: np.ndarray,
    *,
    near_zero_expression_max: float,
) -> dict[str, Any]:
    """Summarize the finite log-expression values across selected samples."""
    array = np.asarray(values)
    if array.ndim != 2 or array.size == 0:
        raise ValueError("Expression values must be a nonempty two-dimensional matrix")
    finite_mask = np.isfinite(array)
    finite = array[finite_mask]
    missing_count = int(np.isnan(array).sum())
    infinite_count = int(np.isinf(array).sum())
    total_count = int(array.size)
    if finite.size:
        zero_count = int(np.sum(finite == 0))
        near_zero_count = int(np.sum((finite >= 0) & (finite <= near_zero_expression_max)))
        finite_count = int(finite.size)
        stats = {
            "minimum_expression": float(np.min(finite)),
            "maximum_expression": float(np.max(finite)),
            "median_expression": float(np.median(finite)),
            "mean_expression": float(np.mean(finite, dtype=np.float64)),
            "n_zero_expression_values": zero_count,
            "zero_expression_fraction": zero_count / finite_count,
            "n_near_zero_expression_values": near_zero_count,
            "near_zero_expression_fraction": near_zero_count / finite_count,
        }
    else:
        stats = {
            "minimum_expression": math.nan,
            "maximum_expression": math.nan,
            "median_expression": math.nan,
            "mean_expression": math.nan,
            "n_zero_expression_values": 0,
            "zero_expression_fraction": math.nan,
            "n_near_zero_expression_values": 0,
            "near_zero_expression_fraction": math.nan,
        }
        finite_count = 0
    return {
        "n_genes": int(array.shape[0]),
        "n_samples": int(array.shape[1]),
        "n_expression_values": total_count,
        "n_finite_values": finite_count,
        "n_missing_values": missing_count,
        "n_infinite_values": infinite_count,
        "missing_fraction": missing_count / total_count,
        "infinite_fraction": infinite_count / total_count,
        "nonfinite_fraction": (missing_count + infinite_count) / total_count,
        **stats,
    }


def tukey_fence_flags(
    values: Sequence[float] | np.ndarray,
    *,
    multiplier: float = 1.5,
) -> tuple[np.ndarray, float, float]:
    """Return strict Tukey-fence flags and fences for a finite metric vector."""
    array = np.asarray(values, dtype=np.float64)
    if array.ndim != 1 or array.size == 0 or not np.isfinite(array).all():
        raise ValueError("Tukey input must be a nonempty vector of finite values")
    if multiplier < 0:
        raise ValueError("Tukey IQR multiplier must be nonnegative")
    q1, q3 = np.quantile(array, [0.25, 0.75], method="linear")
    iqr = float(q3 - q1)
    lower = float(q1 - multiplier * iqr)
    upper = float(q3 + multiplier * iqr)
    flags = (array < lower) | (array > upper)
    return flags, lower, upper


def calculate_sample_metrics(
    values: np.ndarray,
    sample_ids: Sequence[str],
    *,
    selected_expression_threshold: float,
    near_zero_expression_max: float,
    maximum_missing_fraction: float,
    exclude_any_infinite_values: bool,
    exclude_no_finite_values: bool,
    tukey_iqr_multiplier: float,
) -> tuple[list[dict[str, Any]], float, float]:
    """Compute label-free per-sample metrics, deterministic exclusions, and flags."""
    array = np.asarray(values, dtype=np.float64)
    if array.ndim != 2 or array.shape[1] != len(sample_ids):
        raise ValueError("Expression matrix columns must match sample_ids")
    if array.shape[0] == 0 or array.shape[1] == 0:
        raise ValueError("Expression matrix must contain genes and samples")
    if not 0 <= maximum_missing_fraction <= 1:
        raise ValueError("maximum_missing_fraction must be between 0 and 1")

    rows: list[dict[str, Any]] = []
    eligible_indices: list[int] = []
    for column_index, sample_id in enumerate(sample_ids):
        column = array[:, column_index]
        missing_count = int(np.isnan(column).sum())
        infinite_count = int(np.isinf(column).sum())
        finite = column[np.isfinite(column)]
        n_genes = int(column.size)
        missing_fraction = missing_count / n_genes
        infinite_fraction = infinite_count / n_genes
        finite_count = int(finite.size)
        if finite_count:
            median = float(np.median(finite))
            q1, q3 = np.quantile(finite, [0.25, 0.75], method="linear")
            mad = float(np.median(np.abs(finite - median)))
            mean = float(np.mean(finite, dtype=np.float64))
            aggregate_signal = float(np.sum(finite, dtype=np.float64))
            minimum = float(np.min(finite))
            maximum = float(np.max(finite))
            zero_count = int(np.sum(finite == 0))
            near_zero_count = int(np.sum((finite >= 0) & (finite <= near_zero_expression_max)))
            counts = {
                threshold: int(np.sum(finite > threshold)) for threshold in SENSITIVITY_THRESHOLDS
            }
            selected_count = int(np.sum(finite > selected_expression_threshold))
        else:
            median = q1 = q3 = mad = mean = aggregate_signal = minimum = maximum = math.nan
            zero_count = near_zero_count = 0
            counts = dict.fromkeys(SENSITIVITY_THRESHOLDS, 0)
            selected_count = 0
        if exclude_no_finite_values and finite_count == 0:
            included = False
            exclusion_reason = "excluded_no_finite_expression_values"
        elif exclude_any_infinite_values and infinite_count > 0:
            included = False
            exclusion_reason = "excluded_infinite_expression_values"
        elif missing_fraction > maximum_missing_fraction:
            included = False
            exclusion_reason = "excluded_excess_missingness"
        else:
            included = True
            exclusion_reason = "retained_for_qc"
            eligible_indices.append(column_index)
        rows.append(
            {
                "sample_id": sample_id,
                "n_genes_total": n_genes,
                "n_finite_values": finite_count,
                "n_missing_values": missing_count,
                "n_infinite_values": infinite_count,
                "missing_fraction": missing_fraction,
                "infinite_fraction": infinite_fraction,
                "nonfinite_fraction": (missing_count + infinite_count) / n_genes,
                "minimum_expression": minimum,
                "maximum_expression": maximum,
                "mean_expression": mean,
                "median_expression": median,
                "q1_expression": float(q1),
                "q3_expression": float(q3),
                "iqr_expression": float(q3 - q1),
                "median_absolute_deviation": mad,
                "aggregate_expression_signal": aggregate_signal,
                "n_genes_gt_0": counts[0.0],
                "n_genes_gt_1": counts[1.0],
                "n_genes_gt_2": counts[2.0],
                "n_genes_above_selected_threshold": selected_count,
                "n_zero_expression": zero_count,
                "zero_expression_fraction": zero_count / finite_count if finite_count else math.nan,
                "n_near_zero_expression": near_zero_count,
                "near_zero_expression_fraction": (
                    near_zero_count / finite_count if finite_count else math.nan
                ),
                "tukey_outlier_aggregate_signal": None,
                "aggregate_signal_lower_fence": math.nan,
                "aggregate_signal_upper_fence": math.nan,
                "sample_included_in_qc_cohort": included,
                "sample_exclusion_reason": exclusion_reason,
            }
        )

    if not eligible_indices:
        raise ValueError("Sample QC rules excluded every cohort sample")
    eligible_signals = [rows[index]["aggregate_expression_signal"] for index in eligible_indices]
    flags, lower_fence, upper_fence = tukey_fence_flags(
        eligible_signals,
        multiplier=tukey_iqr_multiplier,
    )
    for index, flag in zip(eligible_indices, flags, strict=True):
        rows[index]["tukey_outlier_aggregate_signal"] = bool(flag)
        rows[index]["aggregate_signal_lower_fence"] = lower_fence
        rows[index]["aggregate_signal_upper_fence"] = upper_fence
    return rows, lower_fence, upper_fence


def build_gene_filter_summary(
    gene_ids: Sequence[str],
    values: np.ndarray,
    sample_inclusion_mask: Sequence[bool],
    *,
    sensitivity_thresholds: Sequence[float],
    selected_expression_threshold: float,
    minimum_sample_fraction: float,
) -> GeneFilterResult:
    """Summarize genes and filter by expression prevalence, never by labels."""
    matrix = np.asarray(values, dtype=np.float32)
    included = np.asarray(sample_inclusion_mask, dtype=bool)
    if matrix.ndim != 2 or matrix.shape[0] != len(gene_ids):
        raise ValueError("Gene identifiers must match expression matrix rows")
    if included.ndim != 1 or included.size != matrix.shape[1] or not included.any():
        raise ValueError("sample_inclusion_mask must select at least one matrix column")
    if not 0 < minimum_sample_fraction <= 1:
        raise ValueError("minimum_sample_fraction must be in (0, 1]")
    threshold_set = {float(value) for value in sensitivity_thresholds}
    missing_thresholds = set(SENSITIVITY_THRESHOLDS) - threshold_set
    if missing_thresholds:
        raise ValueError(
            f"Sensitivity thresholds must include {SENSITIVITY_THRESHOLDS}; "
            f"missing {sorted(missing_thresholds)}"
        )

    eligible_values = matrix[:, included]
    finite_mask = np.isfinite(eligible_values)
    finite_count = finite_mask.sum(axis=1).astype(np.int64)
    missing_count = np.isnan(eligible_values).sum(axis=1).astype(np.int64)
    infinite_count = np.isinf(eligible_values).sum(axis=1).astype(np.int64)
    n_samples = int(eligible_values.shape[1])
    required = max(1, math.ceil(minimum_sample_fraction * n_samples))

    counts_by_threshold: dict[float, np.ndarray] = {}
    for threshold in sensitivity_thresholds:
        threshold = float(threshold)
        counts_by_threshold[threshold] = np.sum(
            finite_mask & (eligible_values > threshold), axis=1
        ).astype(np.int64)
    selected_counts = np.sum(
        finite_mask & (eligible_values > selected_expression_threshold), axis=1
    ).astype(np.int64)
    keep_mask = selected_counts >= required
    sensitivity_counts = {
        threshold: int(np.sum(counts_by_threshold[threshold] >= required))
        for threshold in SENSITIVITY_THRESHOLDS
    }

    minimum = np.full(len(gene_ids), math.nan, dtype=np.float64)
    maximum = np.full(len(gene_ids), math.nan, dtype=np.float64)
    means = np.full(len(gene_ids), math.nan, dtype=np.float64)
    medians = np.full(len(gene_ids), math.nan, dtype=np.float64)
    for index in range(len(gene_ids)):
        finite_values = eligible_values[index, np.isfinite(eligible_values[index])]
        if finite_values.size:
            minimum[index] = float(np.min(finite_values))
            maximum[index] = float(np.max(finite_values))
            means[index] = float(np.mean(finite_values, dtype=np.float64))
            medians[index] = float(np.median(finite_values))

    def count(threshold: float) -> np.ndarray:
        return counts_by_threshold[threshold]

    rows: list[dict[str, Any]] = []
    for index, gene_id in enumerate(gene_ids):
        rows.append(
            {
                "gene_id": gene_id,
                "n_samples_in_qc_cohort": n_samples,
                "n_finite_values": int(finite_count[index]),
                "n_missing_values": int(missing_count[index]),
                "n_infinite_values": int(infinite_count[index]),
                "minimum_expression": float(minimum[index]),
                "maximum_expression": float(maximum[index]),
                "mean_expression": float(means[index]),
                "median_expression": float(medians[index]),
                "n_samples_gt_0": int(count(0.0)[index]),
                "fraction_samples_gt_0": float(count(0.0)[index] / n_samples),
                "n_samples_gt_1": int(count(1.0)[index]),
                "fraction_samples_gt_1": float(count(1.0)[index] / n_samples),
                "n_samples_gt_2": int(count(2.0)[index]),
                "fraction_samples_gt_2": float(count(2.0)[index] / n_samples),
                "selected_filter_expression_gt": float(selected_expression_threshold),
                "n_samples_above_selected_threshold": int(selected_counts[index]),
                "fraction_samples_above_selected_threshold": float(
                    selected_counts[index] / n_samples
                ),
                "required_samples_for_selected_filter": required,
                "retained_by_prevalence_filter": bool(keep_mask[index]),
            }
        )
    return GeneFilterResult(
        rows=rows,
        sensitivity_counts=sensitivity_counts,
        required_sample_count=required,
        keep_mask=keep_mask,
    )


def run_expression_qc(
    expression: ExpressionMatrix,
    cohort_records: Sequence[Mapping[str, str]],
    *,
    near_zero_expression_max: float,
    sensitivity_thresholds: Sequence[float],
    selected_expression_threshold: float,
    minimum_sample_fraction: float,
    maximum_missing_fraction: float,
    exclude_any_infinite_values: bool,
    exclude_no_finite_values: bool,
    tukey_iqr_multiplier: float,
) -> QCRunResult:
    """Compute sample/gene QC and independently attach cohort labels for reports."""
    sample_ids = [str(row["sample_id"]) for row in cohort_records]
    if sample_ids != list(expression.sample_ids):
        raise ValueError("Expression sample order must exactly match the classification cohort")
    sample_metrics, lower_fence, upper_fence = calculate_sample_metrics(
        expression.values,
        expression.sample_ids,
        selected_expression_threshold=selected_expression_threshold,
        near_zero_expression_max=near_zero_expression_max,
        maximum_missing_fraction=maximum_missing_fraction,
        exclude_any_infinite_values=exclude_any_infinite_values,
        exclude_no_finite_values=exclude_no_finite_values,
        tukey_iqr_multiplier=tukey_iqr_multiplier,
    )
    for metric_row, cohort_row in zip(sample_metrics, cohort_records, strict=True):
        for field in COHORT_REQUIRED_FIELDS:
            if field != "sample_id":
                metric_row[field] = str(cohort_row[field])

    inclusion_mask = [row["sample_included_in_qc_cohort"] for row in sample_metrics]
    gene_filter = build_gene_filter_summary(
        expression.gene_ids,
        expression.values,
        inclusion_mask,
        sensitivity_thresholds=sensitivity_thresholds,
        selected_expression_threshold=selected_expression_threshold,
        minimum_sample_fraction=minimum_sample_fraction,
    )
    global_stats = calculate_global_expression_stats(
        expression.values,
        near_zero_expression_max=near_zero_expression_max,
    )
    return QCRunResult(
        expression=expression,
        sample_metrics=sample_metrics,
        gene_filter=gene_filter,
        global_stats=global_stats,
        outlier_lower_fence=lower_fence,
        outlier_upper_fence=upper_fence,
    )


def write_tsv(path: Path, rows: Sequence[Mapping[str, Any]], fields: Sequence[str]) -> None:
    """Write a stable TSV schema; nonfinite floating statistics are blank."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(fields), delimiter="\t", extrasaction="ignore"
        )
        writer.writeheader()
        for row in rows:
            output: dict[str, Any] = {}
            for field in fields:
                value = row.get(field, "")
                if value is None:
                    output[field] = ""
                elif isinstance(value, (float, np.floating)) and not math.isfinite(float(value)):
                    output[field] = ""
                elif isinstance(value, (bool, np.bool_)):
                    output[field] = str(bool(value)).lower()
                elif isinstance(value, np.integer):
                    output[field] = int(value)
                elif isinstance(value, np.floating):
                    output[field] = format(float(value), ".10g")
                elif isinstance(value, float):
                    output[field] = format(value, ".10g")
                else:
                    output[field] = value
            writer.writerow(output)
