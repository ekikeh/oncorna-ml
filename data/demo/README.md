# Phase 1G demo data

This folder contains a small, source-derived TCGA-BRCA preview—not synthetic data. It is intended to let a visitor inspect the expression-table schema, sample/patient mappings, normalized PAM50 labels, and provenance without downloading the full cohort or matrix.

- `demo_expression.tsv`: 100 gene rows by 10 sample columns, plus the leading `gene_id` column.
- `demo_metadata.tsv`: one row per expression sample, with sample ID, patient ID, original and normalized PAM50 labels, and the frozen split partition.
- `demo_manifest.json`: deterministic selection rules, source checksums, dimensions, value semantics, and output checksums.

Two samples are included for each of the five available normalized PAM50 classes. Every demo sample comes from the frozen training partition; no validation or test sample is included. This tiny, balanced preview is not population-representative and must not be used for performance estimation or biological conclusions.

The expression values are copied verbatim from the Phase 1E processed matrix and remain `log2(normalized_count + 1)`. No extra transformation, scaling, imputation, or modeling-based feature selection is applied. The full cohort, matrix, and local split manifest are not included here.

For exact selection rules, reproduction prerequisites, limitations and attribution,
see [preview documentation](../../docs/demo_dataset.md) and the
[verified redistribution evidence](../../docs/data_provenance.md#verified-redistribution-evidence-reviewed-2026-10-27).
The expression metadata permits redistribution; NCI and Xena statements support reuse
of the documented open-TCGA preview with attribution. Preserve the upstream versions,
checksum references and selection rules. MIT covers project code only. The sample and
participant identifiers may be linkable and are not anonymous; do not attempt to
reidentify participants. Full-repository snapshots include this preview; installable
demo archives and the minimal image exclude it.
