# OncoRNA-ML

[![CI](https://github.com/ekikeh/oncorna-ml/actions/workflows/ci.yml/badge.svg)](https://github.com/ekikeh/oncorna-ml/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

**A reproducible, leakage-aware study of RNA-seq features for breast-cancer PAM50 subtype classification.**

> **Research/education only.** This project is not a medical device and must not be used for diagnosis, treatment selection, or patient-care decisions.

## Status

**Phase 0 — project setup.** This repository is a starter scaffold. No cohort has been downloaded, no model has been trained, and no performance or biomarker claims are made yet. The first runnable analysis will use a small, deterministic synthetic/demo dataset; full-cohort downloads will remain optional and outside Git.

## Goal

Build an end-to-end, reproducible learning project using public, de-identified breast-cancer data. The primary planned task is multiclass prediction of PAM50 intrinsic subtype. Differential expression and pathway analysis are descriptive/research analyses; a survival analysis is a secondary extension only if endpoints can be harmonized reliably.

The project will prioritize sound validation over a headline score: patient-level splits, feature selection inside cross-validation, macro-F1 and per-class metrics, a locked test set, uncertainty reporting, and explicit limitations.

## Quick start (Phase 0)

Recommended: install **Miniforge** (which provides Conda/Mamba), then run these commands from the repository root:

```bash
conda env create -f environment.yml
conda activate oncorna-ml
python -m pip install -e ".[dev]"
python -m pytest
ruff check .
```

If you use Mamba, replace `conda env create` with `mamba env create`. The environment targets Python 3.11. Windows users should use WSL2 for the smoothest bioinformatics tooling experience.

For this Linux workspace, the same environment can be created with `uv` instead:

```bash
uv python install 3.11
uv venv --python 3.11 .venv
uv pip install --python .venv/bin/python -e ".[dev]"
.venv/bin/python -m pytest
.venv/bin/ruff check .
```

The tests currently check only the project scaffold. The data-processing and modeling pipeline will be added in the next phases.

## Planned data and analysis

- Planned cohort: TCGA-BRCA public, de-identified expression and clinical data; exact source, release, download date, label source, license/terms, and checksums will be recorded before analysis.
- Full cohort files and derived intermediate matrices must **not** be committed. Only a small synthetic or clearly redistributable demo fixture belongs under `data/demo/`.
- Differential expression requires raw gene-level counts. Do not run DESeq2/PyDESeq2 on TPM/FPKM values. A normalized/log-expression source may be used for a separately specified ML analysis, but the source and transformation must be explicit.
- A patient—not an aliquot or sample—is the unit for splitting. Any feature selection or tuning must happen inside the training folds. PAM50 prediction is a reconstruction task, not a claim of discovering new PAM50 biomarkers.

See [`docs/project_charter.md`](docs/project_charter.md), [`docs/data_provenance.md`](docs/data_provenance.md), and [`docs/model_card.md`](docs/model_card.md).

## Repository map

```text
configs/                 analysis settings and random seeds
 data/demo/               small synthetic / redistributable smoke-test fixture only
 data/raw/, data/processed/  local data locations; contents are git-ignored
 docs/                    project charter, provenance, model card
 src/oncorna/              installable Python package
 tests/                    automated checks
 workflow/                 planned Snakemake workflow (added with the first analysis)
 .github/workflows/        lightweight pull-request CI
```

## Public-repository checklist

Repository and citation metadata are initialized for the GitHub account `ekikeh`. Never commit credentials, identifiable information, controlled-access data, or large cohort files. Confirm the data-use and redistribution terms separately from the code license. Keep the disclaimer and limitations visible.

## License and citation

The **code** is licensed under the MIT License. Data and third-party resources retain their own terms; see `docs/data_provenance.md`. Citation metadata is in [`CITATION.cff`](CITATION.cff).
