# Public release scope

## Included now

The tracked repository includes source code, configuration, safe scientific documentation, aggregate development and validation reports, this aggregate final-test summary, source checksums, tests, and a small attributed real-data training preview. The MIT license covers code. The preview and upstream TCGA/Xena data retain their own terms; see [provenance](data_provenance.md) and [demo attribution](demo_dataset.md).

## Excluded

Do not include private fitted-model binaries, full patient-level expression or metadata arrays, frozen local split contents, patient-level validation or test predictions, test authorizations, attempt claims, or sealed audit artifacts in a public package, container, release attachment, or CI input. The final evaluation is terminal and closed. The public aggregate results in the [Phase 3D report](../reports/phase_3d_aggregate_results.md) are the publication boundary.

A public user can inspect the code and reports and run safe tests. The repository does not supply the controlled inputs, fitted model, or authorization needed to reproduce the 169-patient held-out evaluation. The existing 100-gene, 10-sample preview is real training data for structure inspection; it cannot serve as a final-model input or performance demonstration.

## Standalone synthetic demo and package controls

The [standalone synthetic demo](synthetic_demo.md) generates its entire fixture in memory, fits a fresh synthetic classifier, and prints only labeled aggregate results. Its strict ordered-gene validation component accepts an explicit schema and does not load a model. No artifact-loading interface or production inference is provided; caller-supplied models and model-path options are not accepted.

The [packaging controls](packaging.md) provide platform hash locks, an exact wheel/source-archive allowlist and forbidden-content audits, plus a three-module Docker context. The existing real-training preview is excluded from all new archives and images. The repository and installable source archive have different inventories: the full scientific history remains in the repository. A synthetic environment does not reproduce the frozen result or establish clinical use. Inspect actual built contents before publication; ignore rules alone are not proof of exclusion.
