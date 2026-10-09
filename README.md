# OncoRNA-ML

[![CI](https://github.com/ekikeh/oncorna-ml/actions/workflows/ci.yml/badge.svg)](https://github.com/ekikeh/oncorna-ml/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

**A reproducible, leakage-aware study of RNA-seq features for breast-cancer PAM50 subtype classification.**

> **Research/education only.** This project is not a medical device and must not be used for diagnosis, treatment selection, or patient-care decisions.

## Status

**Data-preparation Phases 1A–1G, training-only Phase 2B cross-validation, and one frozen Phase 2C validation evaluation are implemented.** The workflow acquires public UCSC Xena TCGA-BRCA inputs, defines and quality-checks the approved cohort, creates a processed expression matrix, freezes a patient-level train/validation/test split, and commits a compact demo subset. Phase 2B compares a majority baseline with L2 multinomial logistic regression across four training-only folds; Phase 2C fits the selected model on the 506 training patients and evaluates the frozen 169-patient validation set once. See the [Phase 2B report](reports/phase_2b_training_cv.md) and [Phase 2C report](reports/phase_2c_validation.md). The full cohort, processed matrix, frozen split, fitted model, and patient-level predictions remain local-only. No test evaluation, train-plus-validation refit, differential expression, GSEA, or biological claims have been made.

## Goal

Build an end-to-end, reproducible learning project using public, de-identified breast-cancer data. The primary planned task is multiclass prediction of PAM50 intrinsic subtype. Differential expression and pathway analysis are descriptive/research analyses; a survival analysis is a secondary extension only if endpoints can be harmonized reliably.

The project prioritizes sound validation over a headline score: patient-level splits, feature selection inside cross-validation, macro-F1 and per-class metrics, a locked test set, uncertainty reporting, and explicit limitations.

## Quick start (setup and Phase 1A)

The project targets **Python 3.11**. Recommended for a local machine: install Miniforge, then run from the repository root:

```bash
conda env create -f environment.yml
conda activate oncorna-ml
python -m pip install -e ".[dev]"
python -m pytest
ruff check .
python scripts/download_data.py
```

If you use Mamba, replace `conda env create` with `mamba env create`. Windows users should use WSL2 for the smoothest bioinformatics tooling experience.

For this Linux workspace, the same environment can be created with `uv`:

```bash
uv python install 3.11
uv venv --python 3.11 .venv
uv pip install --python .venv/bin/python -e ".[dev]"
.venv/bin/python -m pytest
.venv/bin/ruff check .
.venv/bin/python scripts/download_data.py
```

The download retrieves about 63 MiB of source files, writes SHA-256 checksums to `data/MANIFEST.sha256`, and skips existing files on later runs. Add `--force` to replace them. Tests use small mocked responses; they do not fetch the cohort. Full raw files stay local and are not committed.

The tracked demo tables in [`data/demo/`](data/demo/) can be inspected without downloading the full matrix. Regenerating them requires the local approved cohort, processed matrix, retained-gene list, and frozen `split_v1`; see [`docs/demo_dataset.md`](docs/demo_dataset.md) for the reproduction command and prerequisites.

## Selected data and scientific guardrails

- Cohort: TCGA-BRCA from the public UCSC Xena TCGA hub.
- Expression: `HiSeqV2.gz`, a gene-level, processed legacy TCGA RNA-seq matrix.
- Clinical/sample metadata and PAM50 label: `BRCA_clinicalMatrix`, including `PAM50Call_RNAseq` where available.
- The Xena expression values are log2-transformed normalized expression, **not raw integer counts**. Do not use them directly as DESeq2/PyDESeq2 counts. If a later phase includes count-based differential expression, select and document a suitable raw-count source first.
- Full raw cohort files, the full derived expression matrix, and the local split manifest must **not** be committed. `data/demo/` contains only a small, attributed source-derived preview; it is not a synthetic dataset or a population-representative sample.
- Demo samples come only from the training partition so validation and test examples stay out of the public preview. The fixed gene subset is for compact table inspection, not biological feature selection.
- A patient—not an aliquot or sample—is the unit for the frozen split. Future feature selection and tuning must stay inside training folds. PAM50 prediction is label reconstruction, not a claim of discovering new PAM50 biomarkers.

See [`docs/data_provenance.md`](docs/data_provenance.md) for source URLs, checksums, and usage notes; [`docs/data_splits.md`](docs/data_splits.md) for the frozen split; [`docs/demo_dataset.md`](docs/demo_dataset.md) for the demo rules; and the [`project charter`](docs/project_charter.md) and [`model card`](docs/model_card.md) for scope and limitations.

For the training-only experiment, use Python 3.11 and run `python scripts/run_training_cv.py` after preparing the local `preprocessing_v1` artifacts. The versioned [modeling configuration](configs/modeling_v1.yaml) fixes four stratified folds, the gene filter, scaler, and logistic `C` grid. The command writes local-only results to `data/processed/modeling_v1/` and refuses to overwrite them.

## Repository map

```text
configs/                   project settings and separate Phase 1G demo settings
scripts/                   data acquisition, split, and demo-build commands
data/MANIFEST.sha256       source checksums; no full cohort data
data/demo/                 compact demo expression, metadata, and provenance
data/raw/, data/processed/ local data locations; full files are git-ignored
docs/                      project charter, provenance, split, and demo docs
src/oncorna/               installable Python package
tests/                     automated checks
workflow/                  planned Snakemake workflow (later phase)
.github/workflows/         pull-request CI
```

## License and citation

The **code** is licensed under the MIT License. Data and third-party resources retain their own terms; the Xena source pages do not display a file-level license, and this repository does not claim that the MIT license covers the demo data. The small subset is attributed in [`docs/demo_dataset.md`](docs/demo_dataset.md); see [`docs/data_provenance.md`](docs/data_provenance.md) before reusing it. Citation metadata is in [`CITATION.cff`](CITATION.cff).
