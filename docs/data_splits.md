# Phase 1F — deterministic patient-level data splits

## Purpose and leakage boundary

The approved cohort contains one selected expression sample for each of **844 unique patients**. Splitting by `patient_id` makes the unit of allocation explicit and prevents samples from the same patient from appearing in more than one partition. The split manifest records each patient's selected `sample_id` so the expression sample mapping is preserved.

The assignment uses the normalized PAM50 label only for stratification. It does **not** read expression values, select genes, or use model or performance results. The expression matrix is checked **header-only** to confirm that its sample columns exactly match the approved cohort order. No train/validation/test modeling, differential expression, GSEA, feature selection, survival analysis, or test-set evaluation is performed in Phase 1F.

## Deterministic rounding and stratification

The configured seed is `20261006` (`random_seed` in `configs/default.yaml`). Target counts use the largest-remainder method: floor each `n × proportion`, then allocate leftover patients in descending fractional-remainder order. Ties use the configured fixed precedence **validation, test, train**. For 844 patients and 60/20/20 proportions, this yields **506 train, 169 validation, and 169 test** patients (and one sample per patient).

Stratified assignment uses two fixed stages with scikit-learn's `train_test_split`:

1. Allocate 169 test patients from all 844, stratifying on `pam50_normalized_label`, with random state **20261006**.
2. Allocate 169 validation patients from the remaining 675, stratifying again, with random state **20261007** (the configured seed plus offset 1). The remaining 506 patients form the training partition.

After membership is assigned, patient/sample pairs are serialized in original cohort row order for readability; this does not change membership. The manifest records the random states, scikit-learn version, cohort digest, matrix-header check, counts, and validation results. Re-running the command with identical inputs verifies the existing split and leaves its original timestamp and assignments untouched. If the cohort, seed, configuration, software version, or assignments change, the script refuses to overwrite `split_v1`; use a new version rather than reshuffling the frozen split.

## Actual partition sizes and PAM50 counts

| Normalized PAM50 class | Cohort | Train | Validation | Test |
|---|---:|---:|---:|---:|
| Luminal A | 421 | 253 | 84 | 84 |
| Luminal B | 192 | 115 | 38 | 39 |
| Basal-like | 141 | 85 | 28 | 28 |
| HER2-enriched | 67 | 40 | 14 | 13 |
| Normal-like | 23 | 13 | 5 | 5 |
| **Total** | **844** | **506** | **169** | **169** |

All five normalized PAM50 classes are represented in each partition. Integer allocation means the class fractions are approximate rather than identical across splits.

## Small-class limitation and locked test set

The **Normal-like** class has only 23 patients: 13 are in training and only 5 in each of validation and test. Stratification preserves representation, but cannot remove the substantial uncertainty caused by such a small class. A single test-set Normal-like case is 20% of that class's test count, so class-specific test estimates will be unstable. This limitation should be considered in human review; do not address it by moving patients after viewing model results.

The test partition is **locked during model development**. Do not inspect test-set performance or use test patients for model, feature, or hyperparameter selection. Any eventual evaluation must be a one-time final assessment after development decisions are complete; Phase 1F performs no such evaluation.

## Output and regeneration

`data/processed/split_v1.json` contains the schema version, timestamp, seed and proportions, ordered patient/sample IDs and explicit patient-to-sample mapping for each partition, PAM50 counts, cohort SHA-256 digest, split method, matrix-header validation, and partition validation results. It is Git-ignored with the processed data artifacts.

Create or verify the frozen split from the repository root with:

```bash
python scripts/create_data_splits.py
```

A same-input rerun verifies `split_v1.json` without rewriting it. The script will refuse to replace a frozen split if its provenance or assignment differs.
