# OncoRNA-ML

[![CI](https://github.com/ekikeh/oncorna-ml/actions/workflows/ci.yml/badge.svg)](https://github.com/ekikeh/oncorna-ml/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

**A reproducible, leakage-aware study of RNA-seq features for breast-cancer PAM50 subtype classification.**

> **Research/education only.** This project is not a medical device and must not be used for diagnosis, treatment selection, or patient-care decisions.

## Status

**Phase 1A — reproducible data acquisition.** The Python downloader for public UCSC Xena TCGA-BRCA expression and clinical/PAM50 files is implemented. The files have been retrieved locally and checksummed, but remain ignored by Git. No preprocessing, QC, modeling, differential expression, GSEA, or machine learning has started; no performance or biomarker claims are made. A small synthetic demo will be added in a later phase.

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

## Selected data and scientific guardrails

- Cohort: TCGA-BRCA from the public UCSC Xena TCGA hub.
- Expression: `HiSeqV2.gz`, a gene-level, processed legacy TCGA RNA-seq matrix.
- Clinical/sample metadata and PAM50 label: `BRCA_clinicalMatrix`, including `PAM50Call_RNAseq` where available.
- The Xena expression values are log2-transformed normalized expression, **not raw integer counts**. Do not use them directly as DESeq2/PyDESeq2 counts. If a later phase includes count-based differential expression, select and document a suitable raw-count source first.
- Raw cohort files and derived patient-level matrices must **not** be committed. Only a small synthetic or clearly redistributable demo fixture belongs under `data/demo/`.
- A patient—not an aliquot or sample—is the unit for future splits. Feature selection and tuning must stay inside training folds. PAM50 prediction is label reconstruction, not a claim of discovering new PAM50 biomarkers.

See [`docs/data_provenance.md`](docs/data_provenance.md) for source URLs, date, checksums, and usage notes; also see the [`project charter`](docs/project_charter.md) and [`model card`](docs/model_card.md).

## Repository map

```text
configs/                   analysis settings and random seeds
scripts/download_data.py   public TCGA-BRCA downloader and checksum updater
data/MANIFEST.sha256       checksums only; no cohort data
data/demo/                 future synthetic / redistributable smoke-test fixture
data/raw/, data/processed/ local data locations; full files are git-ignored
docs/                      project charter, provenance, and model card
src/oncorna/                installable Python package
tests/                      automated checks
workflow/                   planned Snakemake workflow (later phase)
.github/workflows/          pull-request CI
```

## License and citation

The **code** is licensed under the MIT License. Data and third-party resources retain their own terms; see [`docs/data_provenance.md`](docs/data_provenance.md). Citation metadata is in [`CITATION.cff`](CITATION.cff).
