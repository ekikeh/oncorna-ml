"""Run Phase 1D expression and sample QC from the repository root.

    python scripts/run_qc.py

Only samples listed in the frozen Phase 1C classification cohort are loaded.
"""

from __future__ import annotations

import argparse
import math
import sys
from collections import Counter
from datetime import datetime
from pathlib import Path
from zoneinfo import ZoneInfo

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_PATH = PROJECT_ROOT / "src"
if str(SRC_PATH) not in sys.path:
    sys.path.insert(0, str(SRC_PATH))

import numpy as np  # noqa: E402
import yaml  # noqa: E402

from oncorna.cohort import PAM50_CLASS_ORDER  # noqa: E402
from oncorna.qc import (  # noqa: E402
    GENE_SUMMARY_FIELDS,
    SAMPLE_METRIC_FIELDS,
    QCRunResult,
    read_classification_cohort,
    read_expression_matrix,
    run_expression_qc,
    write_tsv,
)

DEFAULT_CONFIG = PROJECT_ROOT / "configs" / "default.yaml"
SAMPLE_REASON_LABELS = {
    "retained_for_qc": "Retained for QC",
    "excluded_no_finite_expression_values": "Excluded: no finite expression values",
    "excluded_infinite_expression_values": "Excluded: infinite expression value(s)",
    "excluded_excess_missingness": "Excluded: missingness above configured maximum",
}
FIGURE_NAMES = (
    "sample_expression_distribution.png",
    "expressed_genes_per_sample.png",
    "sample_qc_metric_distributions.png",
    "aggregate_expression_signal_outliers.png",
)


def _load_config(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        config = yaml.safe_load(handle)
    if not isinstance(config, dict) or not isinstance(config.get("qc"), dict):
        raise ValueError(f"Configuration file has no qc section: {path}")
    return config


def _number(value: float, digits: int = 4) -> str:
    if not math.isfinite(float(value)):
        return "NA"
    return f"{float(value):,.{digits}f}"


def _as_percent(value: float) -> str:
    if not math.isfinite(float(value)):
        return "NA"
    return f"{100 * float(value):.4f}%"


def _create_qc_figures(
    result: QCRunResult,
    figure_dir: Path,
    *,
    selected_expression_threshold: float,
) -> list[Path]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure_dir.mkdir(parents=True, exist_ok=True)
    sample_rows = result.sample_metrics
    medians = np.asarray([row["median_expression"] for row in sample_rows], dtype=float)
    q1 = np.asarray([row["q1_expression"] for row in sample_rows], dtype=float)
    q3 = np.asarray([row["q3_expression"] for row in sample_rows], dtype=float)
    included = np.asarray([row["sample_included_in_qc_cohort"] for row in sample_rows], dtype=bool)
    signals = np.asarray([row["aggregate_expression_signal"] for row in sample_rows], dtype=float)
    expressed_counts = np.asarray(
        [row["n_genes_above_selected_threshold"] for row in sample_rows], dtype=float
    )
    outlier_flags = np.asarray(
        [bool(row["tukey_outlier_aggregate_signal"]) for row in sample_rows], dtype=bool
    )

    paths: list[Path] = []

    # Per-sample quantile ribbons show the distribution of all pre-filter genes.
    order = np.argsort(medians, kind="stable")
    ranks = np.arange(1, len(order) + 1)
    fig, ax = plt.subplots(figsize=(11, 5.5))
    ax.fill_between(ranks, q1[order], q3[order], color="#9ecae1", alpha=0.65, label="Q1–Q3")
    ax.plot(ranks, medians[order], color="#08519c", linewidth=1.4, label="Median")
    ax.set(
        xlabel="Samples ordered by median expression",
        ylabel="Expression (log2(normalized_count + 1))",
        title="Per-sample expression distribution before gene filtering",
    )
    ax.legend(frameon=False)
    ax.grid(axis="y", alpha=0.2)
    fig.tight_layout()
    path = figure_dir / FIGURE_NAMES[0]
    fig.savefig(path, dpi=160)
    plt.close(fig)
    paths.append(path)

    fig, ax = plt.subplots(figsize=(9, 5.5))
    ax.hist(expressed_counts[included], bins=30, color="#4c78a8", edgecolor="white")
    median_count = float(np.median(expressed_counts[included]))
    ax.axvline(median_count, color="#e45756", linestyle="--", linewidth=1.5)
    ax.set(
        xlabel=f"Genes with expression > {selected_expression_threshold:g}",
        ylabel="Samples",
        title="Expressed genes per QC-eligible sample",
    )
    ax.text(
        0.98,
        0.95,
        f"Median: {median_count:,.0f}",
        ha="right",
        va="top",
        transform=ax.transAxes,
    )
    ax.grid(axis="y", alpha=0.2)
    fig.tight_layout()
    path = figure_dir / FIGURE_NAMES[1]
    fig.savefig(path, dpi=160)
    plt.close(fig)
    paths.append(path)

    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    qc_panels = (
        (signals[included], "Aggregate expression signal", "Sum of finite log-expression values"),
        (medians[included], "Sample median expression", "log2(normalized_count + 1)"),
        (
            np.asarray([row["zero_expression_fraction"] for row in sample_rows])[included],
            "Zero-expression fraction",
            "Fraction of finite genes equal to zero",
        ),
        (
            np.asarray([row["missing_fraction"] for row in sample_rows])[included],
            "Sample missingness",
            "Fraction of all genes missing",
        ),
    )
    for ax, (values, title, xlabel) in zip(axes.flat, qc_panels, strict=True):
        ax.hist(values, bins=30, color="#72b7b2", edgecolor="white")
        ax.set(title=title, xlabel=xlabel, ylabel="Samples")
        if title == "Sample missingness":
            maximum = float(np.max(values)) if len(values) else 0.0
            ax.set_xlim(0, max(0.01, maximum * 1.1))
            if maximum == 0:
                ax.text(
                    0.98,
                    0.95,
                    "All samples: 0 missing",
                    ha="right",
                    va="top",
                    transform=ax.transAxes,
                )
        ax.grid(axis="y", alpha=0.2)
    fig.suptitle("Sample-level QC metric distributions", y=1.01)
    fig.tight_layout()
    path = figure_dir / FIGURE_NAMES[2]
    fig.savefig(path, dpi=160, bbox_inches="tight")
    plt.close(fig)
    paths.append(path)

    order = np.argsort(signals, kind="stable")
    ranks = np.arange(1, len(order) + 1)
    ordered_flags = outlier_flags[order] & included[order]
    fig, ax = plt.subplots(figsize=(11, 5.5))
    ax.scatter(
        ranks[~ordered_flags],
        signals[order][~ordered_flags],
        s=16,
        color="#7f8c8d",
        alpha=0.65,
        label="Not flagged",
    )
    ax.scatter(
        ranks[ordered_flags],
        signals[order][ordered_flags],
        s=36,
        color="#e45756",
        label="Tukey-fence flag",
        zorder=3,
    )
    ax.axhline(result.outlier_lower_fence, color="#e45756", linestyle="--", linewidth=1)
    ax.axhline(result.outlier_upper_fence, color="#e45756", linestyle="--", linewidth=1)
    flagged_rank = 0
    for rank, sample_index in zip(ranks, order, strict=True):
        if ordered_flags[rank - 1]:
            vertical_offset = 8 if flagged_rank % 2 == 0 else -12
            ax.annotate(
                sample_rows[sample_index]["sample_id"],
                (rank, signals[sample_index]),
                xytext=(7, vertical_offset),
                textcoords="offset points",
                fontsize=8,
            )
            flagged_rank += 1
    ax.set(
        xlabel="Samples ordered by aggregate expression signal",
        ylabel="Aggregate expression signal (sum of log-expression values)",
        title="Exploratory Tukey-fence flags (not automatic exclusions)",
    )
    ax.legend(frameon=False)
    ax.grid(axis="y", alpha=0.2)
    fig.tight_layout()
    path = figure_dir / FIGURE_NAMES[3]
    fig.savefig(path, dpi=160)
    plt.close(fig)
    paths.append(path)
    return paths


def _render_report(
    result: QCRunResult,
    cohort_records: list[dict[str, str]],
    *,
    expression_scale: str,
    gene_filter_expression_gt: float,
    minimum_sample_fraction: float,
    near_zero_expression_max: float,
    maximum_missing_fraction: float,
    exclude_any_infinite_values: bool,
    exclude_no_finite_values: bool,
    tukey_iqr_multiplier: float,
    figure_paths: list[Path],
) -> str:
    stats = result.global_stats
    sample_rows = result.sample_metrics
    gene_filter = result.gene_filter
    excluded_counts = Counter(row["sample_exclusion_reason"] for row in sample_rows)
    total_samples = len(sample_rows)
    final_samples = sum(bool(row["sample_included_in_qc_cohort"]) for row in sample_rows)
    excluded_samples = total_samples - final_samples
    outliers = [
        row
        for row in sample_rows
        if row["sample_included_in_qc_cohort"] and row["tukey_outlier_aggregate_signal"]
    ]
    retained_class_counts = Counter(
        row["pam50_normalized_label"] for row in sample_rows if row["sample_included_in_qc_cohort"]
    )
    all_gene_count = len(result.expression.gene_ids)
    selected_gene_count = int(gene_filter.keep_mask.sum())
    filtered_gene_count = all_gene_count - selected_gene_count
    required = gene_filter.required_sample_count
    threshold_names = {0.0: ">0", 1.0: ">1", 2.0: ">2"}
    date_text = datetime.now(ZoneInfo("Asia/Singapore")).date().isoformat()
    source_cols = result.expression.source_sample_count
    noncohort_cols = result.expression.noncohort_source_sample_count
    ratio = selected_gene_count / all_gene_count if all_gene_count else math.nan

    lines = [
        "# Phase 1D — expression and sample QC",
        "",
        f"**Generated:** {date_text} (Asia/Singapore)",
        "",
        "## Scope and data boundary",
        "",
        f"The input expression unit is `{expression_scale}`. This is a processed, "
        "log-transformed normalized-expression matrix, **not raw integer counts**. No values "
        "are inverse-transformed or described as raw counts. Count-based filtering thresholds "
        "cannot be interpreted directly on this scale; no raw-count threshold was invented.",
        "",
        "Only the sample IDs in `data/processed/classification_cohort.tsv` were loaded. The "
        f"source matrix has **{source_cols:,}** sample columns; **{noncohort_cols:,}** columns "
        "outside the frozen cohort were not read into the QC matrix. No sample outside the "
        "Phase 1C cohort was used.",
        "",
        "No train/validation/test split, model fitting, PAM50-driven gene filtering, "
        "differential expression, GSEA, survival analysis, or ML feature selection was done. "
        "Gene prevalence filtering below uses expression values and sample counts only; PAM50 "
        "labels are attached afterward for cohort counts and outlier reporting.",
        "",
        "## Starting cohort and expression-matrix audit",
        "",
        f"- Starting classification cohort: **{total_samples:,} samples / "
        f"{len({row['patient_id'] for row in cohort_records}):,} unique patients**.",
        f"- Source matrix genes before filtering: **{all_gene_count:,}**.",
        f"- Cohort samples before sample QC: **{total_samples:,}**; source matrix sample columns: "
        f"**{source_cols:,}**.",
        f"- Expression range among finite values: **{_number(stats['minimum_expression'])}** "
        f"to **{_number(stats['maximum_expression'])}**.",
        f"- Median expression across finite gene-by-sample values: "
        f"**{_number(stats['median_expression'])}**.",
        f"- Missing values: **{stats['n_missing_values']:,}** "
        f"({_as_percent(stats['missing_fraction'])}); infinite values: "
        f"**{stats['n_infinite_values']:,}** ({_as_percent(stats['infinite_fraction'])}).",
        f"- Exact zero values: **{stats['n_zero_expression_values']:,}** "
        f"({_as_percent(stats['zero_expression_fraction'])} of finite values).",
        f"- Near-zero values (0 ≤ expression ≤ {near_zero_expression_max:g}): "
        f"**{stats['n_near_zero_expression_values']:,}** "
        f"({_as_percent(stats['near_zero_expression_fraction'])} of finite values).",
        "",
        "The zero/near-zero summaries are descriptive only; they are not additional sample or "
        "gene exclusions.",
        "",
        "## Gene-prevalence filter",
        "",
        "The prevalence denominator is the number of samples remaining after the deterministic "
        "sample-completeness checks below. A gene passes a threshold when its finite expression "
        f"is strictly above that value in at least **{minimum_sample_fraction:.0%}** of eligible "
        f"samples (at least **{required:,} of {final_samples:,}**). Missing and infinite gene "
        "values do not count as expressed.",
        "",
        (
            "| Expression threshold | Genes meeting threshold in ≥"
            f"{minimum_sample_fraction:.0%} of QC samples |"
        ),
        "|---|---:|",
    ]
    for threshold in (0.0, 1.0, 2.0):
        lines.append(
            f"| `{threshold_names[threshold]}` | {gene_filter.sensitivity_counts[threshold]:,} |"
        )
    lines.extend(
        [
            "",
            f"**Selected rule:** expression **> {gene_filter_expression_gt:g}** in at least "
            f"**{minimum_sample_fraction:.0%}** of QC-eligible samples. It retains "
            f"**{selected_gene_count:,} of {all_gene_count:,} genes** "
            f"({ratio:.2%}) and filters {filtered_gene_count:,}. The threshold is a modest "
            "global expression-prevalence screen intended to remove genes absent or nearly absent "
            "in most cohort samples while avoiding a highly restrictive filter. It was fixed "
            "without consulting PAM50 labels and is configurable in `configs/default.yaml`.",
            "",
            "These values remain on the log2-normalized scale. For example, `expression > 1` is "
            "an expression-prevalence criterion on `log2(normalized_count + 1)`; it is not an "
            "integer raw-count threshold and must not be supplied to DESeq2/PyDESeq2.",
            "",
            "## Sample-level QC and exclusions",
            "",
            "Per-sample metrics include aggregate expression signal, the number of genes above "
            "0/1/2 and the selected threshold, mean/median/quartiles/IQR/MAD, zero and near-zero "
            "fractions, missing and infinite values, and minimum/maximum expression.",
            "",
            "`aggregate_expression_signal` is the sum of finite log2-normalized expression values "
            "across the **pre-filter** gene rows for a sample. It is an expression-level summary, "
            "not a read count or library size; it is not a substitute for sequencing depth.",
            "",
            f"Sample exclusion rule: exclude a sample with no finite expression values "
            f"({str(exclude_no_finite_values).lower()}), any infinite value "
            f"({str(exclude_any_infinite_values).lower()}), or missingness greater than "
            f"{maximum_missing_fraction:.1%}. These deterministic checks do not use PAM50. "
            "Tukey outlier flags are exploratory and do **not** automatically exclude samples.",
            "",
            "| Sample disposition | Count |",
            "|---|---:|",
        ]
    )
    for reason in (
        "retained_for_qc",
        "excluded_no_finite_expression_values",
        "excluded_infinite_expression_values",
        "excluded_excess_missingness",
    ):
        lines.append(
            f"| {SAMPLE_REASON_LABELS[reason]} (`{reason}`) | {excluded_counts.get(reason, 0):,} |"
        )
    lines.extend(
        [
            "",
            f"- Samples removed for the listed sample-QC rules: **{excluded_samples:,}**.",
            f"- Final QC cohort: **{final_samples:,} samples / "
            f"{len(retained_class_counts):,} observed PAM50 classes**.",
            "",
            "Class counts are reported for cohort accounting only; they were not used to choose "
            "genes or remove samples:",
            "",
            "| PAM50 class | Final QC samples |",
            "|---|---:|",
        ]
    )
    for label in PAM50_CLASS_ORDER:
        lines.append(f"| `{label}` | {retained_class_counts.get(label, 0):,} |")

    lines.extend(
        [
            "",
            "### Exploratory aggregate-expression outliers",
            "",
            f"Tukey fences use **{tukey_iqr_multiplier:g} × IQR** on "
            "`aggregate_expression_signal` among samples that pass the completeness rule. "
            f"Fences: **{result.outlier_lower_fence:,.3f}** to "
            f"**{result.outlier_upper_fence:,.3f}**. **{len(outliers):,} samples** are flagged; "
            "the flag is not an exclusion.",
            "",
            (
                "| Sample ID | Patient ID | PAM50 | Aggregate expression signal | "
                "Median expression | Genes > selected threshold | Zero fraction | Missingness |"
            ),
            "|---|---|---|---:|---:|---:|---:|---:|",
        ]
    )
    if outliers:
        for row in outliers:
            lines.append(
                f"| `{row['sample_id']}` | `{row['patient_id']}` | "
                f"`{row['pam50_normalized_label']}` | "
                f"{row['aggregate_expression_signal']:,.3f} | "
                f"{row['median_expression']:.4f} | "
                f"{row['n_genes_above_selected_threshold']:,} | "
                f"{_as_percent(row['zero_expression_fraction'])} | "
                f"{_as_percent(row['missing_fraction'])} |"
            )
    else:
        lines.append("| None | — | — | — | — | — | — | — |")

    figures_rel = [path.relative_to(PROJECT_ROOT).as_posix() for path in figure_paths]
    lines.extend(
        [
            "",
            f"All {result.expression.values.shape[0]:,} expression values in each flagged sample "
            "are finite, with no missing or infinite values. The available log2-normalized "
            "matrix does not include "
            "raw read depth, sequencing QC, or independent technical covariates, so it cannot "
            "establish that these low-signal extremes are technical failures rather than "
            "biological variation. **They are retained pending human review; no outlier was "
            "automatically deleted.**",
            "",
            "Human review should compare the flagged samples with source-level sequencing and "
            "clinical/collection metadata before any exclusion decision. Any later exclusion must "
            "be justified by independent technical evidence and encoded as a deterministic rule.",
            "",
            "## Figures",
            "",
        ]
    )
    for path in figures_rel:
        lines.append(f"- [`{path}`]({path.replace('reports/', '../reports/')})")
    lines.extend(
        [
            "",
            "No PCA was run. The figures are descriptive QC visualizations only; none are used "
            "for classification or feature selection.",
            "",
            "## Outputs",
            "",
            "- `data/processed/qc_sample_metrics.tsv`: one row per starting cohort sample, with "
            "cohort identifiers, subtype for reporting, per-sample metrics, outlier flag, and "
            "sample disposition.",
            "- `data/processed/qc_gene_summary.tsv`: one row per source gene, with expression "
            "summaries, prevalence counts/fractions, and the selected gene-filter flag.",
            "- QC figures listed above under `reports/figures/qc/`.",
            "",
            "The TSVs are covered by the repository's `data/processed/` ignore rule. No "
            "processed expression matrix is written or committed.",
            "",
            "## Limitations and scope boundary",
            "",
            "The HiSeqV2 matrix stores `log2(normalized_count + 1)`, not raw integer counts. "
            "Neither the expression-prevalence filter nor the aggregate signal can recover "
            "sequencing depth or support raw-count statistical assumptions. A separate, "
            "documented raw-count source will be required before any DESeq2/PyDESeq2 analysis. "
            "This phase performs no differential expression, GSEA, survival analysis, train/test "
            "split, model training, or ML feature selection.",
            "",
            "Re-run from the repository root with `python scripts/run_qc.py`. The script reads "
            "thresholds and sample-QC choices from `configs/default.yaml`; the QC TSVs are "
            "metadata/statistics outputs only, and no processed expression matrix is copied.",
            "",
        ]
    )
    return "\n".join(lines)


def run(config_path: Path = DEFAULT_CONFIG) -> QCRunResult:
    config = _load_config(config_path)
    qc_config = config["qc"]
    expected_count = int(qc_config["expected_cohort_samples"])
    cohort_records = read_classification_cohort(
        PROJECT_ROOT / qc_config["cohort_path"],
        expected_sample_count=expected_count,
    )
    sample_ids = [row["sample_id"] for row in cohort_records]
    expression = read_expression_matrix(
        PROJECT_ROOT / qc_config["expression_path"],
        sample_ids,
    )
    gene_filter_config = qc_config["gene_filter"]
    sample_filter_config = qc_config["sample_filter"]
    outlier_config = qc_config["outlier_detection"]
    selected_threshold = float(gene_filter_config["expression_gt"])
    minimum_sample_fraction = float(gene_filter_config["minimum_sample_fraction"])
    near_zero_max = float(qc_config["near_zero_expression_max"])
    maximum_missing_fraction = float(sample_filter_config["maximum_missing_fraction"])
    exclude_any_infinite = bool(sample_filter_config["exclude_any_infinite_values"])
    exclude_no_finite = bool(sample_filter_config["exclude_no_finite_values"])
    tukey_multiplier = float(outlier_config["tukey_iqr_multiplier"])
    sensitivity_thresholds = tuple(
        float(value) for value in qc_config["expression_prevalence_thresholds"]
    )

    result = run_expression_qc(
        expression,
        cohort_records,
        near_zero_expression_max=near_zero_max,
        sensitivity_thresholds=sensitivity_thresholds,
        selected_expression_threshold=selected_threshold,
        minimum_sample_fraction=minimum_sample_fraction,
        maximum_missing_fraction=maximum_missing_fraction,
        exclude_any_infinite_values=exclude_any_infinite,
        exclude_no_finite_values=exclude_no_finite,
        tukey_iqr_multiplier=tukey_multiplier,
    )
    sample_metrics_path = PROJECT_ROOT / qc_config["sample_metrics_path"]
    gene_summary_path = PROJECT_ROOT / qc_config["gene_summary_path"]
    report_path = PROJECT_ROOT / qc_config["report_path"]
    figure_dir = PROJECT_ROOT / qc_config["figure_dir"]
    write_tsv(sample_metrics_path, result.sample_metrics, SAMPLE_METRIC_FIELDS)
    write_tsv(gene_summary_path, result.gene_filter.rows, GENE_SUMMARY_FIELDS)
    figure_paths = _create_qc_figures(
        result,
        figure_dir,
        selected_expression_threshold=selected_threshold,
    )
    report = _render_report(
        result,
        cohort_records,
        expression_scale=str(qc_config["expression_scale"]),
        gene_filter_expression_gt=selected_threshold,
        minimum_sample_fraction=minimum_sample_fraction,
        near_zero_expression_max=near_zero_max,
        maximum_missing_fraction=maximum_missing_fraction,
        exclude_any_infinite_values=exclude_any_infinite,
        exclude_no_finite_values=exclude_no_finite,
        tukey_iqr_multiplier=tukey_multiplier,
        figure_paths=figure_paths,
    )
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(report, encoding="utf-8")
    return result


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config",
        type=Path,
        default=DEFAULT_CONFIG,
        help="Project YAML configuration (default: configs/default.yaml)",
    )
    args = parser.parse_args(argv)
    result = run(args.config)
    retained_samples = sum(
        bool(row["sample_included_in_qc_cohort"]) for row in result.sample_metrics
    )
    flagged_samples = sum(
        bool(row["tukey_outlier_aggregate_signal"]) for row in result.sample_metrics
    )
    print(f"Cohort samples before QC: {len(result.expression.sample_ids):,}")
    print(f"Source expression columns: {result.expression.source_sample_count:,}")
    print(
        "Source columns outside the cohort (not loaded): "
        f"{result.expression.noncohort_source_sample_count:,}"
    )
    print(f"Genes before filtering: {len(result.expression.gene_ids):,}")
    print(f"Genes retained by prevalence filter: {int(result.gene_filter.keep_mask.sum()):,}")
    print(f"Samples retained for QC: {retained_samples:,}")
    print(f"Aggregate-expression Tukey flags (not excluded): {flagged_samples:,}")
    print(f"Sample metrics: {PROJECT_ROOT / 'data/processed/qc_sample_metrics.tsv'}")
    print(f"Gene summary: {PROJECT_ROOT / 'data/processed/qc_gene_summary.tsv'}")
    print(f"QC report: {PROJECT_ROOT / 'docs/qc.md'}")
    print(f"QC figures: {PROJECT_ROOT / 'reports/figures/qc'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
