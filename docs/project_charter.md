# Project charter — OncoRNA-ML

**Current status:** Data preparation (Phases 1A–1G), training-only four-fold cross-validation (Phase 2B), frozen validation (Phase 2C), the final 675-patient fit, and the one-time held-out evaluation (Phase 3D) are complete. The [aggregate report](../reports/phase_3d_aggregate_results.md) records approved final statistics. The historical [Phase 2B](../reports/phase_2b_training_cv.md) and [Phase 2C](../reports/phase_2c_validation.md) reports retain their original chronology.

## Research question and completed endpoint

Can a regularized model reconstruct the five existing TCGA-BRCA expression-derived PAM50 calls from gene expression? The completed primary endpoint is multiclass label reconstruction in one research cohort. It does not establish independent biological ground truth, novel biomarkers, external generalization, or clinical validity. The model is not intended for diagnosis, individual prognosis, treatment selection, or triage.

The original charter also proposed biological interpretation, differential expression and pathway enrichment, and optional overall-survival modeling. These remain **uncompleted proposals**, not outputs of the frozen classification experiment. Survival work would require independently reviewed time, event, censoring, and patient-identity harmonization. Count-based differential expression would require appropriate raw counts; the selected Xena expression values are already log-transformed normalized values.

## Eligibility and scientific protocol

The public source, sample inclusion, PAM50 mapping, and one-patient-one-sample resolution are documented in [cohort definition](cohort_definition.md) and [data provenance](data_provenance.md). The approved cohort has 844 primary-tumor patients and the [frozen split](data_splits.md) has 506 training, 169 validation, and 169 test patients. No controlled BAM/FASTQ or attempt to re-identify participants is part of this project.

One shuffled, stratified four-fold split of training patients compared a majority baseline with logistic candidates and selected `C=1`. Fitted gene filtering and scaling stayed within each fitting fold. This was not the nested five-fold/three-fold protocol proposed in the original project outline. The validation estimator fitted 506 patients and retained 17,613 genes. The final estimator fitted the fixed 675 training-plus-validation patients from the ordered 20,530-gene input, retained 17,623 genes, and underwent one terminal 169-patient test assessment. No result-driven retuning follows the test.

The primary headline is macro F1. Final test macro F1 was **0.7934**, balanced accuracy **0.7556**, and accuracy **0.8521** (144/169 correct). Normal-like recall was **2/5**. See the [model card](model_card.md) for class results, data and processing limitations, and research-only use boundaries.

## Reproducibility and publication

The tracked 10-sample by 100-gene [preview](demo_dataset.md) is real, attributed, training-only data for table inspection. It is not a synthetic prediction demo or performance evidence. The historical scientific environment, public implementation, and proposed synthetic package are distinguished in the [reproducibility appendix](reproducibility.md). Full cohort arrays, fitted binaries, patient-level predictions, and protected evaluation records stay outside the proposed public [release scope](release_scope.md). MIT licenses repository code, not upstream data or its preview.

The test is complete and closed. Future interpretation or external validation needs a separate reviewed protocol; it cannot alter the frozen study.
