# Model card — OncoRNA-ML

**Status:** The five-class TCGA-BRCA study is complete through one independently verified, one-time held-out evaluation. This card describes research evidence, not clinical performance. See the [aggregate results report](../reports/phase_3d_aggregate_results.md) and [reproducibility boundaries](reproducibility.md).

## Purpose and intended users

This project is for students and computational oncology researchers studying how gene expression reconstructs existing PAM50 RNA-seq subtype calls. The target is the source `PAM50Call_RNAseq` field, mapped to Luminal A, Luminal B, Basal-like, HER2-enriched, and Normal-like. These labels themselves derive from expression; agreement with them is not independent biological ground truth or new biomarker discovery.

The model must not be used for diagnosis, individual prognosis, treatment selection, clinical triage, or patient care. It is not a medical device and has no clinical or external-cohort validation.

## Data and eligibility

Public UCSC Xena legacy TCGA-BRCA expression and clinical tables supplied the cohort. The [cohort definition](cohort_definition.md) selected one primary tumor sample per patient before requiring a mapped PAM50 call, yielding 844 patients and 844 samples. The [frozen patient split](data_splits.md) contains 506 training, 169 validation, and 169 test patients. The source expression is `log2(normalized_count + 1)`, not raw integer counts; its upstream normalization history cannot be fully reconstructed across partitions. The ordered model input contains 20,530 source genes. The earlier 17,623-gene Phase 1E processed matrix is not the model input.

The [source provenance](data_provenance.md) records URLs, dates, and checksums. MIT covers repository code only. The Xena source pages do not show a file-level license; the [small, source-derived training preview](demo_dataset.md) has its own attribution and reuse caveat. Full patient data and fitted artifacts remain local.

## Model and evaluation

Phase 2B used one shuffled, stratified four-fold split of the 506 training patients. It compared a majority baseline with logistic regression candidates; gene prevalence filtering and standardization were fitted within each fitting fold. Mean fold macro F1 selected `C=1`. This is **not nested cross-validation**: the same four-fold development results selected the hyperparameter and summarize its development performance.

Phase 2C fitted the selected pipeline on 506 training patients, retaining 17,613 genes, and evaluated the frozen 169-patient validation partition once. The final pipeline was then fitted on the fixed 675 training-plus-validation patients, retaining **17,623** genes, and evaluated on the 169-patient held-out test exactly once. Its prevalence filter retains genes with expression strictly greater than 1 in at least 20% of fitting patients, followed by `StandardScaler` and multinomial L2 logistic regression (`C=1`, `lbfgs`, intercept, no class weighting). Learned preprocessing belongs to each fit; the final estimator is distinct from the earlier validation estimator. The test evaluation is terminal and closed. No new tuning follows it.

The original project proposal called for five-fold/three-fold nested CV, calibration, differential expression, pathway analysis, SHAP, and survival modeling. Those analyses were not completed by this study. The actual frozen protocol and its results are reported without rewriting that history.

## Aggregate performance

| Partition and estimator | Macro F1 | Balanced accuracy | Accuracy |
|---|---:|---:|---:|
| Training, four-fold CV (mean of fold macro F1; other columns pooled) | 0.7925 | 0.7891 | 0.8814 |
| Frozen validation, 506-patient fit | 0.7285 | 0.7288 | 0.8284 |
| Held-out test, 675-patient fit | **0.7934** | **0.7556** | **0.8521** |

Training pooled CV macro F1 is **0.7904**, which differs from the mean of four fold scores. The final test has **144 correct classifications among 169 patients**. Test class measures are:

| Class | Precision | Recall | F1 | Support |
|---|---:|---:|---:|---:|
| Luminal A | 0.8556 | 0.9167 | 0.8851 | 84 |
| Luminal B | 0.7297 | 0.6923 | 0.7105 | 39 |
| Basal-like | 1.0000 | 1.0000 | 1.0000 | 28 |
| HER2-enriched | 0.8333 | 0.7692 | 0.8000 | 13 |
| Normal-like | 1.0000 | 0.4000 | 0.5714 | 5 |

The five Normal-like test patients yield only **2/5** correct calls; one patient changes recall by 20 percentage points. Scores across development, validation, and test use different patients and, for validation versus test, different fitted estimators. They should not be read as a performance trend or a confidence interval. No test confidence interval or probability-calibration measure is claimed.

## Limitations and use boundaries

TCGA-BRCA is one research cohort, not a representative clinical population. Expression-derived labels share signal with the predictors. Xena's upstream normalization cannot be independently isolated by this repository, and the prevalence rule originated in earlier cohort-wide QC; therefore an entirely leakage-free end-to-end processing history is not established. Small classes, especially Normal-like, make class estimates unstable. No independent cohort establishes transportability across sites, platforms, or populations. Predicted probabilities have not been calibrated or validated as individual risk estimates.

The public repository provides code, documentation, tests, and a 10-sample by 100-gene **real training-only structural preview**. It does not distribute the fitted model, complete arrays, patient-level predictions, or authorization and seal records. Public readers can inspect the study and run safe tests; they cannot reproduce the closed 169-patient assessment from the public checkout alone. See [release scope](release_scope.md).
