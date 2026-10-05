# Data provenance

This file is a required record, not a claim that data have already been downloaded. Complete it before the first cohort analysis.

## Planned source

- Cohort: TCGA-BRCA (public, de-identified data)
- Expression source and exact assay/file: **to be selected**
- PAM50 label source: **to be selected**
- Clinical source: **to be selected, if used**
- Source URL(s): **to be recorded**
- Release / snapshot date: **to be recorded**
- Download date: **to be recorded**
- Source terms and citation: **to be checked and recorded**
- File checksums (SHA-256): **to be generated after download**

## Data handling rules

- Do not commit full-cohort expression, clinical, or derived patient-level tables. Keep them in ignored local `data/raw/` or `data/processed/` paths.
- Only add a small synthetic or explicitly redistributable fixture under `data/demo/` after documenting its source and license. Never copy real patient IDs into synthetic data.
- Use stable synthetic identifiers such as `SYN-BRCA-0001` for generated demo samples.
- Do not combine counts and normalized-expression sources as if they were interchangeable. Differential expression by DESeq2/PyDESeq2 requires raw counts.
- Do not put API tokens, credentials, or sensitive data in this repository.

## Download manifest

No cohort has been downloaded in Phase 0. When downloads are implemented, record each filename, source URL, retrieval date, byte size, and SHA-256 digest in a generated manifest (for example, `data/MANIFEST.sha256` for redistributable demo files only). Keep the full-cohort manifest local unless its contents are permitted for publication.
