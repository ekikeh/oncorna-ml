# Phase 1E — ML expression matrix preparation

## Scope

This phase materializes the frozen Phase 1C/1D sample-by-gene selection as a reusable, **expression-only** matrix. It does not create data splits, fit models, run differential expression or GSEA, or use PAM50 labels to filter samples or genes. Clinical metadata and subtype labels remain in their existing, separate tables.

## Frozen inputs and selection

The build uses these exact inputs:

- `data/raw/TCGA-BRCA_HiSeqV2.tsv.gz` — source gene-by-sample expression matrix.
- `data/processed/classification_cohort.tsv` — the approved **844-sample** Phase 1C cohort; its row order defines matrix sample-column order. Normal-like samples remain included.
- `data/processed/qc_gene_summary.tsv` — Phase 1D per-gene metrics and prevalence decisions; source row order defines matrix gene-row order.
- `data/processed/qc_sample_metrics.tsv` — Phase 1D sample dispositions and outlier flags; all 844 approved samples passed the sample-completeness rules.
- `configs/default.yaml` — expression units, filter settings, expected dimensions, and output paths.

The expression values are already `log2(normalized_count + 1)`. **The matrix builder copies selected source values verbatim**; it does not apply a second log transform, normalization, centering, scaling, imputation, or other value transformation. The source values are not raw integer counts. The source matrix has 1,218 sample columns; the 374 columns outside the frozen 844-sample cohort are omitted.

The retained genes are exactly the Phase 1D label-independent prevalence-filter result: expression **strictly greater than 1** in at least **20% of 844 samples**, i.e. at least **169 samples**. The expected retained list contains **17,623 of 20,530 source genes**. The builder checks the configured threshold, required sample count, Phase 1D flags, observed source values, row order, and expected dimensions; it fails rather than silently reconciling inconsistent inputs.

The two Tukey-flagged low aggregate-expression samples (`TCGA-C8-A133-01` and `TCGA-D8-A1JS-01`) remain included, as in Phase 1D. The matrix-build manifest records which approved samples carry the outlier flag; that flag is not an exclusion criterion.

## Matrix schema and order

- **Dimensions:** 17,623 gene rows × 844 sample columns (17,623 data rows; 845 TSV columns including the gene-ID column).
- **Orientation:** one gene per row; sample IDs are column headers.
- **First column:** `gene_id`, with source identifiers unchanged.
- **Sample order:** exact row order in `classification_cohort.tsv`, not source-matrix column order.
- **Gene order:** exact source expression-matrix row order, matching `qc_gene_summary.tsv` after applying its Phase 1D retained flag.
- **Values:** verbatim source numeric tokens in `log2(normalized_count + 1)` units.
- **Excluded from the matrix:** PAM50 labels, patient-level clinical fields, QC flags, and other sample metadata.

The companion gene-list TSV makes the selected source gene IDs and prevalence counts explicit. The JSON manifest records SHA-256 digests and paths for all five inputs, the two outputs, units, dimensions, sample and gene order, selected filter, relevant configuration, and the retained outlier flags. It does not embed labels or clinical metadata.

## Outputs and regeneration

Running the command below from the repository root creates:

- `data/processed/ml_expression_matrix.tsv.gz` — full, gzip-compressed TSV expression matrix.
- `data/processed/ml_gene_list.tsv` — selected source gene IDs and Phase 1D prevalence metrics in matrix row order.
- `data/processed/ml_matrix_manifest.json` — machine-readable provenance, configuration, dimensions, order, and file digests.

```bash
python scripts/build_ml_matrix.py
```

Use `python scripts/build_ml_matrix.py --config PATH` to specify another YAML config with the required `qc` and `ml_matrix` settings; keep all three output paths under `data/processed/` so they remain Git-ignored. These **full Phase 1E artifacts are not committed**. Phase 1G separately commits only a compact 100-gene-by-10-training-sample excerpt under `data/demo/`; see [`docs/demo_dataset.md`](demo_dataset.md) for its selection rules, reproduction, and attribution caveat.
