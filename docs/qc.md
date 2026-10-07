# Phase 1D — expression and sample QC

**Generated:** 2026-10-08 (Asia/Singapore)

## Scope and data boundary

The input expression unit is `log2(normalized_count + 1)`. This is a processed, log-transformed normalized-expression matrix, **not raw integer counts**. No values are inverse-transformed or described as raw counts. Count-based filtering thresholds cannot be interpreted directly on this scale; no raw-count threshold was invented.

Only the sample IDs in `data/processed/classification_cohort.tsv` were loaded. The source matrix has **1,218** sample columns; **374** columns outside the frozen cohort were not read into the QC matrix. No sample outside the Phase 1C cohort was used.

No train/validation/test split, model fitting, PAM50-driven gene filtering, differential expression, GSEA, survival analysis, or ML feature selection was done. Gene prevalence filtering below uses expression values and sample counts only; PAM50 labels are attached afterward for cohort counts and outlier reporting.

## Starting cohort and expression-matrix audit

- Starting classification cohort: **844 samples / 844 unique patients**.
- Source matrix genes before filtering: **20,530**.
- Cohort samples before sample QC: **844**; source matrix sample columns: **1,218**.
- Expression range among finite values: **0.0000** to **20.9784**.
- Median expression across finite gene-by-sample values: **7.5845**.
- Missing values: **0** (0.0000%); infinite values: **0** (0.0000%).
- Exact zero values: **2,492,687** (14.3859% of finite values).
- Near-zero values (0 ≤ expression ≤ 0.1): **2,493,131** (14.3884% of finite values).

The zero/near-zero summaries are descriptive only; they are not additional sample or gene exclusions.

## Gene-prevalence filter

The prevalence denominator is the number of samples remaining after the deterministic sample-completeness checks below. A gene passes a threshold when its finite expression is strictly above that value in at least **20%** of eligible samples (at least **169 of 844**). Missing and infinite gene values do not count as expressed.

| Expression threshold | Genes meeting threshold in ≥20% of QC samples |
|---|---:|
| `>0` | 18,615 |
| `>1` | 17,623 |
| `>2` | 16,775 |

**Selected rule:** expression **> 1** in at least **20%** of QC-eligible samples. It retains **17,623 of 20,530 genes** (85.84%) and filters 2,907. The threshold is a modest global expression-prevalence screen intended to remove genes absent or nearly absent in most cohort samples while avoiding a highly restrictive filter. It was fixed without consulting PAM50 labels and is configurable in `configs/default.yaml`.

These values remain on the log2-normalized scale. For example, `expression > 1` is an expression-prevalence criterion on `log2(normalized_count + 1)`; it is not an integer raw-count threshold and must not be supplied to DESeq2/PyDESeq2.

## Sample-level QC and exclusions

Per-sample metrics include aggregate expression signal, the number of genes above 0/1/2 and the selected threshold, mean/median/quartiles/IQR/MAD, zero and near-zero fractions, missing and infinite values, and minimum/maximum expression.

`aggregate_expression_signal` is the sum of finite log2-normalized expression values across the **pre-filter** gene rows for a sample. It is an expression-level summary, not a read count or library size; it is not a substitute for sequencing depth.

Sample exclusion rule: exclude a sample with no finite expression values (true), any infinite value (true), or missingness greater than 5.0%. These deterministic checks do not use PAM50. Tukey outlier flags are exploratory and do **not** automatically exclude samples.

| Sample disposition | Count |
|---|---:|
| Retained for QC (`retained_for_qc`) | 844 |
| Excluded: no finite expression values (`excluded_no_finite_expression_values`) | 0 |
| Excluded: infinite expression value(s) (`excluded_infinite_expression_values`) | 0 |
| Excluded: missingness above configured maximum (`excluded_excess_missingness`) | 0 |

- Samples removed for the listed sample-QC rules: **0**.
- Final QC cohort: **844 samples / 5 observed PAM50 classes**.

Class counts are reported for cohort accounting only; they were not used to choose genes or remove samples:

| PAM50 class | Final QC samples |
|---|---:|
| `Luminal A` | 421 |
| `Luminal B` | 192 |
| `Basal-like` | 141 |
| `HER2-enriched` | 67 |
| `Normal-like` | 23 |

### Exploratory aggregate-expression outliers

Tukey fences use **1.5 × IQR** on `aggregate_expression_signal` among samples that pass the completeness rule. Fences: **123,504.941** to **139,717.690**. **2 samples** are flagged; the flag is not an exclusion.

| Sample ID | Patient ID | PAM50 | Aggregate expression signal | Median expression | Genes > selected threshold | Zero fraction | Missingness |
|---|---|---|---:|---:|---:|---:|---:|
| `TCGA-C8-A133-01` | `TCGA-C8-A133` | `Luminal A` | 122,591.286 | 7.0096 | 15,583 | 19.5616% | 0.0000% |
| `TCGA-D8-A1JS-01` | `TCGA-D8-A1JS` | `Luminal A` | 122,812.177 | 6.9602 | 15,565 | 17.6717% | 0.0000% |

All 20,530 expression values in each flagged sample are finite, with no missing or infinite values. The available log2-normalized matrix does not include raw read depth, sequencing QC, or independent technical covariates, so it cannot establish that these low-signal extremes are technical failures rather than biological variation. **They are retained pending human review; no outlier was automatically deleted.**

Human review should compare the flagged samples with source-level sequencing and clinical/collection metadata before any exclusion decision. Any later exclusion must be justified by independent technical evidence and encoded as a deterministic rule.

## Figures

- [`reports/figures/qc/sample_expression_distribution.png`](../reports/figures/qc/sample_expression_distribution.png)
- [`reports/figures/qc/expressed_genes_per_sample.png`](../reports/figures/qc/expressed_genes_per_sample.png)
- [`reports/figures/qc/sample_qc_metric_distributions.png`](../reports/figures/qc/sample_qc_metric_distributions.png)
- [`reports/figures/qc/aggregate_expression_signal_outliers.png`](../reports/figures/qc/aggregate_expression_signal_outliers.png)

No PCA was run. The figures are descriptive QC visualizations only; none are used for classification or feature selection.

## Outputs

- `data/processed/qc_sample_metrics.tsv`: one row per starting cohort sample, with cohort identifiers, subtype for reporting, per-sample metrics, outlier flag, and sample disposition.
- `data/processed/qc_gene_summary.tsv`: one row per source gene, with expression summaries, prevalence counts/fractions, and the selected gene-filter flag.
- QC figures listed above under `reports/figures/qc/`.

The TSVs are covered by the repository's `data/processed/` ignore rule. No processed expression matrix is written or committed.

## Limitations and scope boundary

The HiSeqV2 matrix stores `log2(normalized_count + 1)`, not raw integer counts. Neither the expression-prevalence filter nor the aggregate signal can recover sequencing depth or support raw-count statistical assumptions. A separate, documented raw-count source will be required before any DESeq2/PyDESeq2 analysis. This phase performs no differential expression, GSEA, survival analysis, train/test split, model training, or ML feature selection.

Re-run from the repository root with `python scripts/run_qc.py`. The script reads thresholds and sample-QC choices from `configs/default.yaml`; the QC TSVs are metadata/statistics outputs only, and no processed expression matrix is copied.
