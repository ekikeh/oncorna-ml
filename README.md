# OncoRNA-ML

[![CI](https://github.com/ekikeh/oncorna-ml/actions/workflows/ci.yml/badge.svg)](https://github.com/ekikeh/oncorna-ml/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

**TCGA-BRCA expression-based reconstruction of five existing PAM50 subtype labels.**

> **Research and education only.** This is not a diagnostic, prognostic, treatment-selection, or patient-care tool.

## Research question and status

Can a regularized model reconstruct existing expression-derived PAM50 calls from TCGA-BRCA gene-expression measurements under a frozen patient-level evaluation? The completed study includes cohort preparation (Phases 1A–1G), four-fold training development CV (Phase 2B), one frozen validation assessment (Phase 2C), a final 675-patient fit, and one independently verified held-out test assessment (Phase 3D). The held-out evaluation is terminal and closed. [Historical reports](reports/README.md) retain what was known at each phase.

The five classes are Luminal A, Luminal B, Basal-like, HER2-enriched, and Normal-like. Their source labels come from `PAM50Call_RNAseq`, itself based on expression. The study measures label reconstruction within one cohort, not independent biological ground truth or clinical utility.

## Data and method

The public UCSC Xena TCGA-BRCA legacy `HiSeqV2` expression matrix and `BRCA_clinicalMatrix` supplied the inputs and labels. Values are `log2(normalized_count + 1)`, **not raw integer counts**. Count-based DESeq2/PyDESeq2 analysis requires a separately documented appropriate raw-count source. [Provenance](docs/data_provenance.md) records source URLs, dates, and checksums.

The [cohort rules](docs/cohort_definition.md) selected one eligible primary tumor sample per patient, yielding 844 patients. The [patient-level frozen split](docs/data_splits.md) was 506 training / 169 validation / 169 test. The approved input schema contained 20,530 genes in source order. The fitted pipeline applied a fitting-population prevalence filter (expression > 1 in at least 20% of patients), standardization, and multinomial L2 logistic regression with `C=1`.

[Training CV](reports/phase_2b_training_cv.md) used one stratified four-fold split to compare a majority baseline and logistic candidates and select `C`; it was **not nested CV**. The [validation](reports/phase_2c_validation.md) estimator fitted 506 patients and retained 17,613 genes. The final estimator fitted the fixed 675 training-plus-validation patients and retained 17,623 genes before the one-time test. These are separate fits.

## Verified aggregate results

| Evaluation | Macro F1 | Balanced accuracy | Accuracy |
|---|---:|---:|---:|
| Training CV (mean fold macro F1; other measures pooled) | 0.7925 | 0.7891 | 0.8814 |
| Frozen validation | 0.7285 | 0.7288 | 0.8284 |
| Final held-out test | **0.7934** | **0.7556** | **0.8521** |

Training pooled CV macro F1 is 0.7904, distinct from the mean of four fold scores. The final test correctly classified **144/169** patients; Normal-like recall was **2/5**. See the [aggregate final report](reports/phase_3d_aggregate_results.md) and [model card](docs/model_card.md) for class results and interpretation. No confidence interval, calibration measure, or external-cohort score is claimed.

## Limitations

TCGA-BRCA is one research cohort. Expression-derived labels share signal with expression predictors. Five Normal-like patients in the test make its recall unstable. The upstream Xena normalization history cannot be reconstructed across partitions, and the gene-prevalence threshold arose in earlier cohort-wide QC. An entirely leakage-free end-to-end history is not established. No independent external-cohort validation, clinical calibration, differential expression, GSEA, SHAP, or survival analysis has been completed.

## Quick start: public code and tests

Use Python 3.11. With Miniforge or Mamba, from the repository root:

```bash
conda env create -f environment.yml
conda activate oncorna-ml
python -m pip install -e ".[dev]"
python -m pytest
ruff check .
```

These checks use repository tests; they do not run the completed 169-patient evaluation. The `pyproject.toml` dependencies are ranges, not a lock of the historical scientific environment. The tracked [CI workflow](.github/workflows/ci.yml) runs lint and tests. [Reproducibility notes](docs/reproducibility.md) distinguish the historical environment, aggregate evidence, and public implementation.

For source acquisition, see [data provenance](docs/data_provenance.md) and `scripts/download_data.py`. Downloading and processing the full cohort is separate from this quick start. The tracked [preview](data/demo/README.md) is **10 real training samples by 100 genes**, selected for table inspection. It is neither synthetic nor representative and cannot be passed to the final fitted model. Its [selection and attribution](docs/demo_dataset.md) are documented.

Run `python -m oncorna.synthetic_demo` for a separate offline [synthetic teaching demonstration](docs/synthetic_demo.md), labeled `synthetic_demo_only`. It generates invented data, trains a fresh demonstration model, and prints aggregate results and a confusion matrix. Exact dependency locks and broader CI/packaging changes remain Phase 4D work. The fitted real model, complete arrays, patient-level predictions, and evaluation authorizations remain local; the [release scope](docs/release_scope.md) explains why the public checkout cannot reproduce the closed held-out assessment.

## Repository map

```text
configs/                   versioned scientific settings
scripts/                   acquisition, cohort, split, modeling commands
data/MANIFEST.sha256       source checksums
data/demo/                 attributed real training-only preview
data/raw/, data/processed/ local data and fitted artifacts; ignored by Git
docs/                      cohort, provenance, model card, reproducibility
reports/                   historical and aggregate scientific reports
src/oncorna/               installable Python package
tests/                     automated checks
workflow/                  future Snakemake plan
.github/workflows/         CI
```

## License, citation, and contributions

The [MIT license](LICENSE) covers project code. TCGA/Xena data and other third-party resources retain their own terms; the Xena source pages do not display a file-level license, and MIT does not relicense the preview. Review [data provenance](docs/data_provenance.md) and [demo attribution](docs/demo_dataset.md) before reuse. Citation metadata is in [CITATION.cff](CITATION.cff). Contributions follow [CONTRIBUTING.md](CONTRIBUTING.md).
