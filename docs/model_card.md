# Model card — OncoRNA-ML

**Status:** Phase 2B training-only cross-validation and one frozen Phase 2C validation evaluation completed. The selected pipeline was fitted on all 506 training patients and evaluated once on 169 validation patients; the test set remains locked. Results and limitations are documented in the [Phase 2B report](../reports/phase_2b_training_cv.md) and [Phase 2C report](../reports/phase_2c_validation.md). They are not clinical performance claims.

## Intended use

Educational/research reproduction of a PAM50 subtype label-reconstruction workflow using public, de-identified breast-cancer data. A future model would be an educational research artifact only.

## Out-of-scope use

Diagnosis, prognosis for an individual, treatment selection, clinical triage, or any patient-care decision. This project is not a medical device and has not been clinically validated.

## Planned evaluation (not yet performed)

- **Target:** normalized PAM50 intrinsic-subtype labels derived from the source `PAM50Call_RNAseq` field and documented in [`docs/cohort_definition.md`](cohort_definition.md).
- **Unit of separation:** patient, not aliquot/sample; the frozen train/validation/test design is documented in [`docs/data_splits.md`](data_splits.md).
- **Validation:** any future development must keep the test partition locked and confine learned preprocessing, feature selection, and tuning to training data.
- **Metrics:** macro-F1 as the primary metric, with class-wise precision/recall, balanced accuracy, confusion matrix, and uncertainty if supported by the data.
- **Baselines:** dummy classifier, followed by regularized multinomial logistic regression before considering more complex models.

The Phase 2B development and Phase 2C validation metrics are reported separately. No test result or deployable model artifact has been produced; the fitted pipeline is local-only.

## Data and limitations

The approved TCGA-BRCA classification cohort contains **844 patients and 844 selected samples**. The normalized label counts are Luminal A (421), Luminal B (192), Basal-like (141), HER2-enriched (67), and Normal-like (23). The Phase 1E expression matrix contains 17,623 genes by 844 samples and stores the source's processed `log2(normalized_count + 1)` values; it is not a raw-count matrix. Provenance, inclusion rules, QC, and source limitations are detailed in the Phase 1 data documentation.

The Phase 1G preview contains only 10 frozen-training samples (two per class) and 100 fixed-position genes. It exists for inspecting table structure and source-value semantics, not for model development or performance estimates. TCGA is a research cohort and is not population-representative. PAM50 labels derived from RNA expression can share signal with expression predictors; any future findings must be framed as label reconstruction, not independent biomarker discovery. The upstream Xena source pages do not show a file-level license; the demo attribution and reuse caveat are documented in [`docs/demo_dataset.md`](demo_dataset.md).

## Update policy

Complete this card only after a properly locked evaluation. Report dataset version, split definition, software version, results with uncertainty, subgroup limitations, known failure modes, and model artifact checksum. Do not publish a clinical performance claim.
