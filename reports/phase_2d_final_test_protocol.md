# Phase 2D — final-test protocol

**Status:** the scientific protocol was approved for infrastructure implementation in Phase 2D-B. This approval does not authorize the real 675-patient fit or access to the 169 test patients' expression or labels. The working base is commit `7bb92028391f6132983b580c09e7672c78d715d9`. The companion configuration is [`configs/final_test_v1.yaml`](../configs/final_test_v1.yaml); its Phase 2D-B byte SHA-256 is `1cadc4bce4d891f904a06409ba12794afa854abf8733276dc4484e2f074baa36`. The executable source revision and configuration checksum require human review before any separate execution authorization.

## 1. Decision and prediction task

The final estimator **will be refitted on the union of the original 506 training and 169 validation patients**, exactly 675 distinct patients and one approved primary-tumor sample per patient. The original frozen test partition remains 169 patients. Membership and patient-to-sample mappings in `split_v1` must remain unchanged. The union preserves its recorded patient/sample order; no new random split, resampling, exclusion, or duplicate handling is permitted. From the published frozen split counts, the 675 fitting labels comprise Luminal A 337, Luminal B 153, Basal-like 113, HER2-enriched 54, and Normal-like 18. These are planning counts from existing documentation, not newly inspected test labels.

The task remains single-label reconstruction of five existing `PAM50Call_RNAseq`-derived classes from the Xena TCGA-BRCA expression measurements. The fixed output order is **Luminal A, Luminal B, Basal-like, HER2-enriched, Normal-like**. Patient/sample IDs align inputs and outputs but never become features. The Phase 2C validation metrics apply to a **different, 506-trained estimator**; the 675-trained final estimator has no separate untouched validation estimate.

## 2. Locked model and preprocessing

| Item | Frozen specification |
|---|---|
| Input | `preprocessing_v1` all-gene view, 20,530 source gene IDs in source order, values already `log2(normalized_count + 1)` |
| Gene screen | Existing `GenePrevalenceFilter`: finite expression **strictly > 1** in at least `ceil(0.20 × 675) = 135` fitting patients; fit afresh on the 675 patients and apply its mask unchanged to test |
| Scaling | `StandardScaler`, fit after filtering on those same 675 patients only |
| Classifier | Multinomial L2 logistic regression; `C=1`, `solver=lbfgs`, intercept enabled, `class_weight=None`, `tol=0.0001`, `max_iter=1000` |
| Excluded steps | No imputer, supervised feature selector, PCA, calibration, class weighting, decision threshold, extra model, or hyperparameter search |

The final retained-gene count and ordered gene-list digest are **unknown until the approved 675-patient fit**. The 17,613-gene 506-training mask from Checkpoint 3B/Phase 2C is a provenance comparison only and must not be supplied as an input mask. The same existing filter implementation and class mapping must be used. No model coefficient or gene-importance interpretation is part of final-test execution.

## 3. One-time execution protocol for a later approved phase

1. **Freeze implementation before access.** Implement the runner and focused tests against synthetic fixtures without opening the real test expression or labels. Review the exact code diff, configuration, test plan, and allowed outputs. Commit or otherwise identify an immutable source revision and record the SHA-256 of `configs/final_test_v1.yaml`. Obtain explicit human authorization for the *one* real test run after review. Do not treat approval of this draft as run authorization.
2. **Preflight without test values.** Confirm Git revision and a clean or fully accounted-for working tree. Check the SHA-256 values listed in the companion config for source expression, cohort, all-gene array/schema, preprocessing manifest, frozen split, Phase 2B/2C manifests, and Phase 2B/2C configuration files. Check the frozen split's original 506/169/169 counts, patient/sample-map integrity, disjointness, completeness, and the exact 675-member train-plus-validation union. Verify ordered gene IDs, units, no duplicate IDs, fitting-label mapping, and finite *fitting* expression. Do not inspect test expression or test labels while selecting or fitting the model. Abort before fitting on any discrepancy; never regenerate the split or overwrite frozen artifacts.
3. **Reserve the run.** Atomically create a new `data/processed/final_test_v1/` directory with exclusive creation and refusal if it already exists. Write a start/provenance record before any fit. A directory left by a failed or partial attempt remains an audit artifact; do not delete it to make a repeat look like the first run.
4. **Fit on 675 only.** Materialize only the exact approved train-plus-validation rows from the original 20,530-gene view. Fit a single complete `GenePrevalenceFilter → StandardScaler → LogisticRegression` pipeline. Verify the filter saw precisely those 675 sample IDs, its prevalence count is 135, the scaler saw 675 rows, class membership/order is correct, and all retained features are finite. Record retained-gene count and ordered digest, solver iterations, fit duration, memory/runtime notes, all warnings, software versions, and model checksum. A failed fit, nonfinite value, feature mismatch, or `ConvergenceWarning` stops **before** test access. Record the known scikit-learn `penalty='l2'` deprecation warning without silently changing the approved estimator.
5. **Access test once after a valid fit.** Materialize only the 169 frozen test sample rows in recorded order, check schema and finiteness, and call the fitted pipeline's prediction and probability methods once on that batch. Align the five probability columns by the classifier's actual `classes_` to the approved class order. Confirm exactly one row per test patient, unique sample IDs, valid predicted labels, finite probabilities in `[0,1]`, and row sums near one. Save the prediction/probability artifact locally. Only then read the fixed test labels for scoring; they cannot affect any learned parameter. Do not search thresholds, recalibrate, refit, or re-predict after seeing scores.
6. **Score and seal.** Compute only the metrics below using the frozen label order. Save complete local metrics, raw and row-normalized confusion matrices, prediction distribution, model, training-gene list, run/failure manifest, and SHA-256 of each output. Verify output hashes and the test prediction table's 169 unique frozen patient IDs. Produce a Git-safe aggregate report after reviewing it for patient-level content. Mark the run complete without overwriting or rerunning it.

The existing 506-trained Phase 2C pipeline must not be evaluated on test as an extra comparator. The single final-test estimator is the newly fitted 675-trained pipeline.

## 4. Fixed metrics and report layout

**Primary:** unweighted macro F1 across all five classes, with the fixed class order and `zero_division=0`. **Supplementary:** balanced accuracy, overall accuracy, per-class precision/recall/F1 and support, a raw confusion matrix, a row-normalized confusion matrix with true classes as rows, predicted class counts, and aligned five-class probability outputs. Probabilities are uncalibrated model scores; no probability threshold or calibration claim is allowed. Report all five classes even if a class receives no predictions.

Machine-readable metrics use JSON numbers at full Python float64 serialization precision. The raw confusion matrix has true classes as rows and predicted classes as columns, both in the fixed class order. Row-normalized machine values are fractions in `[0, 1]`, obtained by dividing each cell by that true class's support. Human-readable reports show metrics to four decimal places and row percentages to one decimal place. Probability columns follow the fixed class order by mapping the estimator's actual `classes_`; each row must sum to one within an absolute tolerance of `1e-9`. Every raw confusion row sum must match true-class support, every column sum must match predicted-class count, and both totals must equal the test population. These choices are fixed in `configs/final_test_v1.yaml` and do not change the model or metrics.

The final aggregate report should contain: (a) data/provenance and fit diagnostics; (b) a headline metric table with exact denominators; (c) a five-row class table; (d) raw and row-normalized five-by-five confusion tables or heatmaps; (e) predicted class distribution; (f) a descriptive comparison with Phase 2B training CV and the Phase 2C validation result; (g) warnings/failures and scientific limitations. Do not rank or choose a new model from those comparisons. The Normal-like test support is only five under the documented frozen split, so one case changes its recall by 20 percentage points. Explicitly discuss Normal-like and the previously weak HER2-enriched validation class. Do not present fold spread as a test confidence interval or infer clinical utility from a point estimate.

## 5. Failure and deviation policy

- **Before test access:** if any checksum, membership, schema, finite-value, fitting, class-map, or convergence gate fails, preserve a failure manifest, stop, and request human review. A corrected runner or environment may be proposed under a new recorded version; no parameter or class change is automatic.
- **After any test access:** preserve partial artifacts and an exact access/failure log. Do not automatically repeat the test, alter the model, or use observed predictions/metrics to decide a fix. Any mechanical recovery requires a documented root cause, a narrowly scoped repair, a new versioned output directory, and explicit human authorization. The original attempt remains part of the report.
- **After a completed score:** no performance-based threshold, feature, model, exclusion, or hyperparameter changes. Unexpectedly low or high performance is a result to report, not a reason to rerun or retune. Any departure from this protocol must be visibly labeled as a deviation; it cannot be described as the prespecified final evaluation.

## 6. Provenance and artifact policy

The locked settings are in `configs/final_test_v1.yaml`. The future execution must record the exact code revision and config file hash, original split hash and partition-identity hashes, 675 fitting patient/sample digest, all-gene schema and source hashes, fitted ordered gene-list digest, model hash, library and Python versions, solver iterations, warning classes/messages, timestamps, output hashes, and execution command. The source data and split hashes are already declared in the config but must be independently checked at run time. The environment should match the documented Python 3.11/scikit-learn behavior or be explicitly reviewed before test access.

The new `data/processed/final_test_v1/` directory will be **local only**. It should hold the fitted pipeline, training manifest/gene list, one patient-level test prediction/probability table, aggregate metrics, confusion matrices, and a provenance/failure manifest. These files, full data, and the frozen split must not be committed. Git may contain the reviewed configuration, runner and synthetic tests, an aggregate Markdown report, and aggregate plots only after checking that no patient identifiers, patient-level probabilities, model binary, credentials, or source-derived full-cohort values are included. The project code license does not confer a data redistribution license.

## 7. Scientific interpretation limits

PAM50 labels are derived from RNA expression, as are the predictors, so successful reconstruction does not establish an independent biomarker or a clinical diagnostic. The Normal-like class is scarce, and Phase 2C HER2-enriched recall was only 6/14. The `>1` in 20% prevalence rule was inherited from earlier cohort-wide QC, even though the final mask will be fitted only on 675 patients. Xena's upstream normalization cannot be reconstructed sufficiently to prove that no information crossed the eventual split before this project received the matrix. TCGA is one selected research cohort, with no external cohort validation or population transportability assessment. The 675-trained estimator is different from the 506-trained validation estimator and has no separate untouched validation set; only the one final test can provide its prespecified internal held-out estimate. None of these results should guide patient care.

## 8. Human approval checklist

Complete **all** items before any real test expression or label access:

- [ ] Confirm the original frozen split, cohort, source, and all-gene hashes match the config; confirm the 506/169/169 memberships and exact 675-member union without changing them.
- [ ] Review and approve the fixed `C=1` pipeline, 135-of-675 prevalence requirement, five-class order, no-extra-search rule, and the metrics above.
- [ ] Implement and review the one-shot runner and synthetic tests; record the exact source revision, config checksum, environment versions, and output schema.
- [ ] Confirm the dedicated final-test output directory does not exist, can be created exclusively, and will preserve partial failures without overwrite.
- [ ] Confirm the runner prevents test expression/labels from entering fitting and delays label access until predictions are fixed.
- [ ] Confirm the failure/deviation policy, patient-level local-only policy, aggregate report format, and no-retuning rule.
- [ ] Obtain **explicit human approval for the single final test run** after the preceding items are reviewable and fixed.

**Review decision:** scientific method approved for implementation in Phase 2D-B. The real final-test execution remains pending separate human approval. The Phase 2D-B runner contains an explicit code-level block on real execution.
