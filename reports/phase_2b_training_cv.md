# Phase 2B — training-only cross-validation report

**Run:** 2026-10-09, 13:27–13:29 Singapore time. **Scope:** development on the 506 frozen training patients only. No model was fitted on all 506 patients; no validation or test expression or performance was evaluated. The outcome is reconstruction of existing expression-derived PAM50 labels, not independent biomarker discovery or clinical validation.

## Protocol and integrity

The run used the 20,530-gene `preprocessing_v1` view, Python 3.11.9, NumPy 2.4.6, pandas 2.3.3, and scikit-learn 1.9.1. The full input hashes matched the frozen split (`33e464aa…`), cohort (`d3b117a3…`), raw expression (`263bf672…`), all-gene array (`692eb3c5…`), and gene schema (`fd37450a…`), as well as the preprocessing manifest. The frozen split and historical matrix provenance checks passed. All 506 training rows were finite, gene IDs were unique and in source order, patient/sample/label alignment passed, and training counts were 253 Luminal A, 115 Luminal B, 85 Basal-like, 40 HER2-enriched, and 13 Normal-like. The two retained low-expression outliers were not excluded.

The local out-of-fold file contains exactly 506 predictions for each of five candidates (2,530 rows). Its patient-ID set equals the frozen training set, with **zero validation and zero test patient-ID overlap**. It is Git-ignored. The runner selects only training rows from the memory-mapped expression array after checking artifact hashes and schema; the historical 17,623-gene matrix and the full-training 17,613-gene list do not define CV features.

One shuffled, stratified four-fold split with seed `20261009` was reused across all candidates. Within each logistic fit, `GenePrevalenceFilter(>1, 20%)` and `StandardScaler` were fitted on the fold's fitting patients only. The classifier used multinomial L2 logistic regression with `lbfgs`, intercept, no class weighting, tolerance `0.0001`, and maximum 1,000 iterations. Each candidate was run sequentially in one process. The dummy used the most-frequent training-fold class. No imputer was needed because the training matrix has no missing values.

| Fold | Fitting patients | Held-out patients | Fitting class counts, in fixed label order | Held-out Luminal A | Luminal B | Basal-like | HER2-enriched | Normal-like | Genes retained in logistic fit |
|---:|---:|---:|---|---:|---:|---:|---:|---:|---:|
| 1 | 379 | 127 | 189 / 87 / 63 / 30 / 10 | 64 | 28 | 22 | 10 | 3 | 17,626 |
| 2 | 379 | 127 | 190 / 86 / 64 / 30 / 9 | 63 | 29 | 21 | 10 | 4 | 17,611 |
| 3 | 380 | 126 | 190 / 86 / 64 / 30 / 10 | 63 | 29 | 21 | 10 | 3 | 17,614 |
| 4 | 380 | 126 | 190 / 86 / 64 / 30 / 10 | 63 | 29 | 21 | 10 | 3 | 17,606 |

Each fold had disjoint fitting and held-out patients. Every training patient was held out exactly once per candidate. The gene counts were independently refitted in every fold and were the same across `C` candidates because filtering is label-independent and the fold membership did not change.

## Development results

**Primary ranking:** arithmetic mean of the four fold macro F1 scores; exact ties go to smaller `C`. Standard deviation is the sample standard deviation of four fold scores. Pooled metrics score the 506 out-of-fold predictions together and are distinct from the arithmetic mean of fold metrics.

| Candidate | Fold macro F1: 1 / 2 / 3 / 4 | Mean ± SD | Pooled macro F1 | Pooled balanced accuracy | Pooled accuracy | Fit time, four folds |
|---|---|---:|---:|---:|---:|---:|
| Majority dummy | .1340 / .1326 / .1333 / .1333 | .1333 ± .0006 | .1333 | .2000 | .5000 | 0.2 s |
| Logistic `C=.001` | .8002 / .7649 / .7659 / .7998 | .7827 ± .0200 | .7819 | .7438 | .8794 | 74.6 s |
| Logistic `C=.01` | .8002 / .7649 / .7659 / .8006 | .7829 ± .0202 | .7822 | .7448 | .8794 | 51.5 s |
| Logistic `C=.1` | .8002 / .6740 / .7653 / .7907 | .7575 ± .0576 | .7583 | .7265 | .8735 | 19.9 s |
| **Logistic `C=1`** | .8881 / .7426 / .7427 / .7968 | **.7925 ± .0687** | **.7904** | **.7891** | **.8814** | 9.1 s |

`C=1` is selected because its mean fold macro F1 (.7925) is highest; no tie-break was needed. Its advantage over `C=.01` is about .0096 in mean fold macro F1, while its fold spread is wider. This is an approved-protocol selection, not evidence of a robust superiority claim. The selected model correctly called 6 of 13 Normal-like patients in pooled out-of-fold predictions; that class remains a major limitation.

The following table gives pooled out-of-fold **precision / recall / F1** for each class. Support is fixed across candidates: Luminal A 253, Luminal B 115, Basal-like 85, HER2-enriched 40, Normal-like 13. Exact unrounded values and per-fold class-wise metrics are in the local `candidate_results.json`.

| Candidate | Luminal A | Luminal B | Basal-like | HER2-enriched | Normal-like |
|---|---|---|---|---|---|
| Dummy | .500 / 1.000 / .667 | 0 / 0 / 0 | 0 / 0 / 0 | 0 / 0 / 0 | 0 / 0 / 0 |
| `C=.001` | .867 / .957 / .910 | .838 / .765 / .800 | 1.000 / .965 / .982 | .829 / .725 / .773 | .800 / .308 / .444 |
| `C=.01` | .870 / .953 / .909 | .832 / .774 / .802 | 1.000 / .965 / .982 | .829 / .725 / .773 | .800 / .308 / .444 |
| `C=.1` | .867 / .953 / .908 | .827 / .748 / .785 | .988 / .976 / .982 | .806 / .725 / .763 | .750 / .231 / .353 |
| **`C=1`** | .914 / .921 / .917 | .852 / .800 / .825 | .923 / .988 / .955 | .775 / .775 / .775 | .500 / .462 / .480 |

Pooled out-of-fold confusion matrices below use rows as true labels and columns as predicted labels. The fixed order for both axes is **Luminal A, Luminal B, Basal-like, HER2-enriched, Normal-like**. The matching per-class support above is the row total; row-normalized percentages can be recovered by dividing each row by its support.

```text
Dummy:    [253,0,0,0,0] [115,0,0,0,0] [85,0,0,0,0] [40,0,0,0,0] [13,0,0,0,0]
C=.001:  [242,9,0,1,1] [26,88,0,1,0] [1,0,82,2,0] [3,8,0,29,0] [7,0,0,2,4]
C=.01:   [241,10,0,1,1] [25,89,0,1,0] [1,0,82,2,0] [3,8,0,29,0] [7,0,0,2,4]
C=.1:    [241,10,0,1,1] [27,86,0,2,0] [0,0,83,2,0] [3,8,0,29,0] [7,0,1,2,3]
C=1:     [233,12,1,3,4] [18,92,1,3,1] [0,0,84,1,0] [1,4,3,31,1] [3,0,2,2,6]
```

For the selected `C=1`, the **row-normalized** confusion matrix below expresses each true class as percentages across predicted classes, in the same fixed order. Rounding may make a row differ slightly from 100%.

| True class | Luminal A | Luminal B | Basal-like | HER2-enriched | Normal-like |
|---|---:|---:|---:|---:|---:|
| Luminal A | 92.1% | 4.7% | 0.4% | 1.2% | 1.6% |
| Luminal B | 15.7% | 80.0% | 0.9% | 2.6% | 0.9% |
| Basal-like | 0% | 0% | 98.8% | 1.2% | 0% |
| HER2-enriched | 2.5% | 10.0% | 7.5% | 77.5% | 2.5% |
| Normal-like | 23.1% | 0% | 15.4% | 15.4% | 46.2% |

There were **no failed fits and no convergence warnings**. All 16 logistic fits emitted a scikit-learn 1.9.1 `FutureWarning` because its explicit `penalty='l2'` argument is deprecated in that version; the warning and iteration count were recorded for each fold. It does not indicate failed convergence. The four logistic candidates took 155.1 seconds of measured fitting time in total. End-to-end elapsed time, including input checks and output writing, was about 160 seconds. These measurements are machine-specific.

## Reproduction and artifacts

From the repository root with the Python 3.11 environment active, run `python scripts/run_training_cv.py`. The command refuses to overwrite `data/processed/modeling_v1/`. To rerun, deliberately configure a new versioned output directory; the same inputs and settings should reproduce fold membership, predictions, and metrics, while timestamps and fit times may change. Synthetic tests also compare repeated executions with identical seed and data. The local output contains:

- `fold_assignments.json` — exact fit/hold patient IDs and indices, SHA-256 `4f41a83c…`;
- `candidate_results.json` — complete fold and pooled metrics, warnings, gene counts, SHA-256 `815202bc…`;
- `oof_predictions.tsv` — training-patient-level predictions, SHA-256 `10b31048…`;
- `manifest.json` — input/config/output hashes, versions, Git revision at run time, class counts, seed, selected `C`, and selection rule.

All four files stay local under Git-ignored `data/processed/modeling_v1/`. No full expression data, split file, patient-level prediction table, or fitted model is committed. The run manifest records base Git commit `d9c87af7d65ad2fb9360181a5c99f90322404ef6` and a dirty tree because Phase 2B source files were present at run time. The committed source revision is recorded in the final handoff, not retroactively inserted into the immutable local manifest.

## Verification, interpretation, and next checkpoint

The focused modeling tests cover patient-disjoint and reproducible folds, complete out-of-fold coverage, fixed label order, held-out preprocessing isolation, hyperparameter tie-break, warning recording, and failure recording. The complete suite passed **72 tests** when pytest's temporary directory was placed inside the writable workspace. Ruff lint passed. Ruff formatting passed for all new files; the repository-wide format check still identifies the pre-existing `tests/test_download_data.py` discrepancy. The first full test attempt used the Windows sandbox's default temporary directory and failed on permission errors; its rerun inside the workspace passed without source changes.

These are **development** metrics and may be optimistic because the same CV results selected `C`. Four folds share training patients, so their standard deviation is descriptive variability, not a confidence interval. Normal-like has only 13 training patients, with 3–4 held out per fold; its class-wise estimates are particularly unstable. PAM50 labels and predictors both derive from RNA expression. The Xena normalization history cannot be fully reconstructed across partitions, and the fixed gene-prevalence rule originated in earlier cohort-wide QC. One selected TCGA cohort provides no external generalizability or clinical validation.

**Ready for Phase 2C human review:** the selected `C=1` and training-only pipeline are documented. Phase 2C must decide whether and when to fit all 506 training patients and open the frozen validation set. This Phase 2B checkpoint did neither and did not touch the frozen test evaluation.
