# Phase 3D — aggregate held-out results

**Status:** The approved one-time held-out evaluation completed and was independently verified. Approval `ONCORNA-PHASE3D-TEST-APPROVAL-001`; attempt `7589bebb-7a2f-4279-bc90-cc44cf72ddda`. The evaluation is terminal and closed. This report records approved aggregate statistics and prior protocol summaries; it does not derive new results from sealed or patient-level artifacts.

## Scientific task and frozen protocol

The task reconstructs five existing expression-derived `PAM50Call_RNAseq` labels from legacy TCGA-BRCA log-transformed normalized expression. The [cohort](../docs/cohort_definition.md) has one eligible primary-tumor sample per patient, 844 patients total. The [frozen patient split](../docs/data_splits.md) assigned 506 to training, 169 to validation, and 169 to test.

[Phase 2B](phase_2b_training_cv.md) used one stratified four-fold training CV design to select `C=1` from a majority baseline and L2 multinomial logistic-regression candidates. The filter and scaler were fitted afresh within each fold. This was not nested CV: the fold scores were also used for selection. [Phase 2C](phase_2c_validation.md) fitted the selected pipeline on 506 training patients, retained 17,613 genes, and evaluated validation once. The final estimator fitted the fixed 675 training-plus-validation patients. Its ordered input schema contained 20,530 genes, its fitted prevalence filter retained 17,623, and it used `StandardScaler` and L2 multinomial logistic regression with `C=1`. That final fit is distinct from the validation fit. The 169-patient test was then evaluated once.

## Aggregate results

| Evaluation | Macro F1 | Balanced accuracy | Accuracy |
|---|---:|---:|---:|
| Training CV: mean of four fold macro F1 scores; pooled other measures | 0.7925 | 0.7891 | 0.8814 |
| Training CV: pooled out-of-fold predictions | 0.7904 | 0.7891 | 0.8814 |
| Frozen validation, 506-patient estimator | 0.7285 | 0.7288 | 0.8284 |
| Final held-out test, 675-patient estimator | **0.7934** | **0.7556** | **0.8521** |

The mean of four macro F1 scores (0.7925) and macro F1 computed from pooled out-of-fold predictions (0.7904) are different summaries. The final held-out estimator correctly classified **144 of 169** patients.

| Test class | Precision | Recall | F1 | Support |
|---|---:|---:|---:|---:|
| Luminal A | 0.8556 | 0.9167 | 0.8851 | 84 |
| Luminal B | 0.7297 | 0.6923 | 0.7105 | 39 |
| Basal-like | 1.0000 | 1.0000 | 1.0000 | 28 |
| HER2-enriched | 0.8333 | 0.7692 | 0.8000 | 13 |
| Normal-like | 1.0000 | 0.4000 | 0.5714 | 5 |

Normal-like recall is **2/5**. A single patient changes this recall by 20 percentage points. These partition scores are descriptive; they come from different patients, and the validation and test estimators have different fitting populations. No confidence intervals, AUROC, AUPRC, calibration statistic, or new test metric is supplied here.

## Interpretation and boundaries

These results measure reconstruction of expression-derived labels within one TCGA-BRCA cohort. They do not establish independent subtype ground truth, external transportability, biomarker discovery, or clinical validity. Source values are `log2(normalized_count + 1)`, not raw integer counts. Upstream normalization cannot be fully reconstructed across partitions, and the prevalence threshold has a historical cohort-wide QC context. An entirely leakage-free end-to-end history is not claimed. Probabilities are uncalibrated and must not be presented as individual risk estimates.

This public aggregate report excludes patient-level examples, predictions, fitted binaries, authorizations, claims, and sealed audit records. See the [model card](../docs/model_card.md), [reproducibility appendix](../docs/reproducibility.md), and [release scope](../docs/release_scope.md).
