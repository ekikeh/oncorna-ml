# Public release scope

## Included now

The tracked repository includes source code, configuration, safe scientific documentation, aggregate development and validation reports, this aggregate final-test summary, source checksums, tests, and a small attributed real-data training preview. The MIT license covers code. The preview and upstream TCGA/Xena data retain their own terms; see [provenance](data_provenance.md) and [demo attribution](demo_dataset.md).

## Excluded

Do not include private fitted-model binaries, full patient-level expression or metadata arrays, frozen local split contents, patient-level validation or test predictions, test authorizations, attempt claims, or sealed audit artifacts in a public package, container, release attachment, or CI input. The final evaluation is terminal and closed. The public aggregate results in the [Phase 3D report](../reports/phase_3d_aggregate_results.md) are the publication boundary.

A public user can inspect the code and reports and run safe tests. The repository does not supply the controlled inputs, fitted model, or authorization needed to reproduce the 169-patient held-out evaluation. The existing 100-gene, 10-sample preview is real training data for structure inspection; it cannot serve as a final-model input or performance demonstration.

## Proposed later package

A separately reviewed release can add an invented synthetic fixture, deterministic demo, strict ordered-gene input checks, and a safe prediction contract that requires an explicit trusted model identity. Its exact dependency pins, CI run, and any Docker image will apply to that synthetic demonstration. They will not imply that the synthetic outputs reproduce the frozen scientific result or establish clinical use. Inspect built archives and container contents before publication; ignore rules alone are not proof of exclusion.
