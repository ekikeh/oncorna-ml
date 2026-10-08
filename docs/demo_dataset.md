# Phase 1G — compact, reproducible demo dataset

## Purpose and scope

The demo is a small, source-derived TCGA-BRCA excerpt for inspecting table shape, expression-value semantics, sample/patient mappings, normalized PAM50 labels, and provenance. It is not synthetic, is not population-representative, and is not evidence of biological or clinical performance. Phase 1G performs no model training or evaluation, differential expression, enrichment analysis, survival analysis, or biological interpretation.

A clean clone can inspect the tracked demo files without downloading the full cohort or expression matrix. Recreating them requires the local Phase 1C/1E/1F inputs listed below.

## Files, dimensions, and alignment

| File | Format | Contents |
|---|---|---|
| `data/demo/demo_expression.tsv` | UTF-8, tab-separated | 100 gene rows; `gene_id` followed by 10 sample-ID columns (101 total columns). |
| `data/demo/demo_metadata.tsv` | UTF-8, tab-separated | 10 sample rows with `sample_id`, `patient_id`, `pam50_original_label`, `pam50_normalized_label`, and `split`. |
| `data/demo/demo_manifest.json` | JSON | Source checksums, selection rules and selected IDs, dimensions, value semantics, and output checksums. It does not contain the full cohort or split. |

Metadata rows are in the same order as the sample columns in the expression table. The `gene_id` values are in the same order as the selected positions in the Phase 1E retained-gene list. Every sample ID, patient ID, original label, and normalized label is copied from the approved classification cohort; the metadata `split` value is `train` for every row.

## Deterministic sample selection

The builder reads the existing, Git-ignored `split_v1.json`; it never generates, reshuffles, or writes the split. It verifies the frozen split's schema, cohort mapping, cohort checksum, seed, and validation status, then uses **training patients only**. This keeps the validation and locked test examples out of the public preview and does not expose or tune on test data.

For each normalized PAM50 class in the canonical order `Luminal A`, `Luminal B`, `Basal-like`, `HER2-enriched`, `Normal-like`:

1. Take the patients assigned to `train` by the frozen split and join them to the approved cohort by the explicit patient/sample mapping.
2. Compute `SHA-256("oncorna-ml-phase-1g-demo-v1|20261006|<patient_id>")` for each candidate.
3. Sort by the hexadecimal digest, breaking a hypothetical tie by `patient_id`, and take the first **two** patients in that class.
4. Write selected records in canonical class order and then hash-rank order.

This uses only the frozen partition, patient IDs, normalized PAM50 class, and the existing seed; it does not inspect expression values. All five classes have at least two training patients, so the demo includes **2 samples per class, 10 total**. A future input that cannot supply two training patients in every class fails clearly rather than silently dropping or substituting a class.

## Deterministic gene selection

The source list is the **17,623-gene Phase 1E retained-gene list**, in its existing matrix/source row order. For `N = 17,623` and `M = 100`, the builder selects zero-based positions

```text
index(j) = floor(((2*j + 1) * N) / (2*M)), for j = 0, 1, ..., M-1
```

and retains those gene identifiers in source order. This transparent positional rule spreads the compact table across the existing list. It does not inspect the expression values, PAM50 labels, model results, outcome associations, or validation/test performance. It is a display-size rule—not biological feature selection. The full 17,623-gene retained list was established in Phase 1D/1E using the previously documented label-independent prevalence filter; Phase 1G does not re-rank or re-filter it.

## Expression-value semantics

`demo_expression.tsv` is copied from `data/processed/ml_expression_matrix.tsv.gz`, the Phase 1E matrix. Values remain **`log2(normalized_count + 1)`** from the selected public Xena HiSeqV2 source; they are not raw integer counts. The builder selects the source string tokens for the chosen gene/sample intersections and writes them unchanged. It performs no additional log transform, normalization, imputation, scaling, rounding, or model-based selection. The Phase 1E matrix and gene-list checksums are verified against their manifest before extraction.

## Reproduction

From the repository root, with Python 3.11 and the project dependencies installed, run:

```bash
python scripts/build_demo_dataset.py
```

The script reads `configs/default.yaml` for the existing Phase 1C/1E/1F paths and seed, and `configs/demo.yaml` for the Phase 1G sample/gene counts and output paths. These are deliberately separate: the frozen Phase 1F split records a checksum of `configs/default.yaml`, so Phase 1G settings must not alter that file or invalidate the frozen split.

Required local files:

- `data/processed/classification_cohort.tsv`
- `data/processed/split_v1.json` (the existing frozen split; do not regenerate it)
- `data/processed/ml_expression_matrix.tsv.gz`
- `data/processed/ml_gene_list.tsv`
- `data/processed/ml_matrix_manifest.json`

The builder validates the matrix and gene-list checksums, sample/gene dimensions, source ordering, and split/cohort mapping. It writes only the three demo outputs under `data/demo/`; identical inputs produce byte-identical outputs. If the required local inputs are unavailable, inspect the committed demo tables rather than generating substitute or fabricated values.

## Limitations and data attribution

- This is a tiny, deliberately balanced structural preview (two samples per class), not a random sample of the cohort's class frequencies and not a representative clinical or population sample.
- Only training-partition samples are included. The demo cannot be used to evaluate a model, tune a threshold, select genes, estimate performance, or make biological or clinical claims.
- The source is the public/open TCGA-BRCA `TCGA.BRCA.sampleMap/HiSeqV2` expression resource hosted by [UCSC Xena](https://xenabrowser.net/datapages/?dataset=TCGA.BRCA.sampleMap%2FHiSeqV2&host=https%3A%2F%2Ftcga.xenahubs.net), with PAM50 labels from the approved Xena/TCGA clinical source. See [`docs/data_provenance.md`](data_provenance.md) for exact inputs and citations.
- Xena's pages do not display a file-level license for these source files. The [GDC FAQ](https://gdc.cancer.gov/about-gdc/gdc-faqs) describes accredited use of open-access GDC data, but does not establish a separate file-specific Xena license or explicitly resolve republication of this excerpt. This repository provides source attribution and makes no claim that its MIT code license applies to the data. Review current upstream terms before broader redistribution; see [`docs/data_provenance.md`](data_provenance.md).
