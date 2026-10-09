# Phase 2C — frozen validation evaluation

**Status:** one approved evaluation completed on 2026-10-09 at 22:06 Singapore time. The five existing PAM50 expression-derived labels were reconstructed from RNA expression. These results are research/education evidence of label reconstruction, not independent biomarker discovery or clinical validation. The test set remains locked.

## Frozen inputs and execution boundary

The checkout began clean at approved commit `f58abd350d0d49295249c29f6fbf898620987481`. The one-shot run used that HEAD and the unchanged, checksummed Phase 2B `modeling_v1` configuration and local candidate results, which selected logistic `C=1`. The Phase 2B output checksums passed. The `preprocessing_v1` manifest, source expression, all-gene array and schema, cohort, historical matrix provenance, and frozen `split_v1` all passed integrity checks. The split remained **506 train / 169 validation / 169 test**. All training and validation expression values were finite, IDs were unique, labels aligned to the frozen patient/sample pairs, and the 20,530 source genes remained in the approved order. No historical prefitted gene list or 17,623-gene matrix was used as a model input.

Only 506 training rows fitted the pipeline. The validation input contained exactly its 169 frozen patients, with **zero train or test patient overlap**; the saved prediction table has one row per validation patient and five probabilities aligned to the approved order. Each probability vector sums to one. The runner indexed no test expression rows, performed no test prediction, and did not train on combined train-plus-validation patients.

The complete pipeline was `GenePrevalenceFilter(expression > 1 in at least 20% of fitting samples) → StandardScaler → multinomial L2 LogisticRegression(C=1, solver=lbfgs, intercept=True, class_weight=None, tol=.0001, max_iter=1000)`. Gene filtering and scaling were fitted on training patients only. The full-training threshold required 102 of 506 patients and retained **17,613 genes**. Its ordered gene-list SHA-256 is `d4ad384eb78ae3336427405c22bff9bb4ab8b9703f01cb87192994b527c88452`. Scaling recorded 506 fitting samples. The solver stopped after **13 iterations**; fitting took **1.96 seconds**. There was **no convergence warning or failed fit**. scikit-learn 1.9.1 emitted one `FutureWarning` for the approved explicit `penalty='l2'` argument, which is deprecated in that version; this did not alter the selected configuration.

## Validation results

The primary metric, **macro F1, was 0.7285**. Balanced accuracy was **0.7288** and overall accuracy was **0.8284**. Metrics below use the fixed class order **Luminal A, Luminal B, Basal-like, HER2-enriched, Normal-like**, with undefined precision handled as zero.

| True subtype | Support | Precision | Recall | F1 | Predicted as this subtype |
|---|---:|---:|---:|---:|---:|
| Luminal A | 84 | .8941 | .9048 | .8994 | 85 |
| Luminal B | 38 | .7714 | .7105 | .7397 | 35 |
| Basal-like | 28 | .8235 | 1.0000 | .9032 | 34 |
| HER2-enriched | 14 | .6000 | .4286 | .5000 | 10 |
| Normal-like | 5 | .6000 | .6000 | .6000 | 5 |

The raw confusion matrix uses true classes as rows and predicted classes as columns, in the fixed order above.

| True class | Luminal A | Luminal B | Basal-like | HER2-enriched | Normal-like |
|---|---:|---:|---:|---:|---:|
| Luminal A | 76 | 4 | 0 | 2 | 2 |
| Luminal B | 8 | 27 | 1 | 2 | 0 |
| Basal-like | 0 | 0 | 28 | 0 | 0 |
| HER2-enriched | 1 | 4 | 3 | 6 | 0 |
| Normal-like | 0 | 0 | 2 | 0 | 3 |

The row-normalized matrix shows the percentage of each true class assigned to each predicted class. Rounding may make row totals differ slightly from 100%.

| True class | Luminal A | Luminal B | Basal-like | HER2-enriched | Normal-like |
|---|---:|---:|---:|---:|---:|
| Luminal A | 90.5% | 4.8% | 0% | 2.4% | 2.4% |
| Luminal B | 21.1% | 71.1% | 2.6% | 5.3% | 0% |
| Basal-like | 0% | 0% | 100% | 0% | 0% |
| HER2-enriched | 7.1% | 28.6% | 21.4% | 42.9% | 0% |
| Normal-like | 0% | 0% | 40.0% | 0% | 60.0% |

Normal-like recall is **3/5**. One patient changes that class's validation recall by 20 percentage points, so its apparent increase from Phase 2B should not be read as a robust improvement. HER2-enriched recall is **6/14**, with several calls assigned to Luminal B or Basal-like. These are descriptive error counts only; they did not trigger a change to the approved model.

## Comparison with Phase 2B

| Measure | Phase 2B selected `C=1`, training CV | Frozen validation | Difference, validation minus CV |
|---|---:|---:|---:|
| Macro F1 | .7925 mean of four folds | .7285 | −.0641 |
| Balanced accuracy | .7891 pooled out of fold | .7288 | −.0604 |
| Accuracy | .8814 pooled out of fold | .8284 | −.0530 |

The macro-F1 comparison uses a **mean of fold scores** on the CV side; the other CV metrics pool one out-of-fold prediction per training patient. These are different estimands and cohorts, so differences are descriptive. Phase 2B CV also selected `C`, which can make its development estimate optimistic. The validation result was evaluated once and was not used to tune, select, or refit any model.

## Reproducibility and review boundary

Python 3.11.9, NumPy 2.4.6, pandas 2.3.3, scikit-learn 1.9.1, and joblib 1.6.0 were recorded. The local run spans 22:06:15–22:06:17 Singapore time after preflight; the measured fit was 1.96 seconds. Its manifest records the approved Git SHA, input/config hashes, selected `C`, training gene digest, software versions, timestamps, and hashes for each output. It also records a dirty working tree at run time because the Phase 2C implementation files were being prepared for this local commit. The saved pipeline, training manifest, 169 patient-level predictions/probabilities, gene list, aggregate metrics, and both confusion matrices reside only in Git-ignored `data/processed/validation_v1/`. The output directory was new and is now protected against overwrite or repeat evaluation. The fitted pipeline SHA-256 is `57644d6be73d2c9a09e0d25251b193787c0f40d7d528bb181577ba73ee98c38d`; the validation prediction table SHA-256 is `cad527dbea25b5657f3b5522fa60bdbf56392df99a0a2717db4038b8f11ba7c2`. The binary model must be loaded only from this trusted local source with a compatible environment.

Focused tests cover training-only fitting, isolation from held-out expression, no test expression indexing, five-class probability alignment, one prediction per validation patient, metric calculation, immutable output behavior, deterministic synthetic results, and Phase 2B provenance verification. The complete Python 3.11 suite passed **77 tests** from a workspace-local temporary directory. Ruff lint and formatting checks passed for the new files. The repository-wide formatter still reports the known unrelated discrepancy in `tests/test_download_data.py`.

PAM50 labels and predictors both derive from RNA expression, the Normal-like class is very small, and this single legacy TCGA cohort does not establish transportability. The upstream Xena normalization process cannot be fully reconstructed across partitions, and the gene-prevalence threshold originated in earlier cohort-wide QC. This validation result does not support clinical use. **The 169-patient test set remains untouched for a later, separately approved final assessment.** No 675-patient fit, extra model search, or coefficient interpretation was performed in Phase 2C.
