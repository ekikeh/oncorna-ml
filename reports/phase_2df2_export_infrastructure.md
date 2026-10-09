# Phase 2D-F2 — fitting-only expression export infrastructure

**Status:** Synthetic implementation and review only, 2026-10-10. No real fitting artifact was created, no real all-gene expression or test labels were opened, no real 675-patient model was fitted, and no final test was executed.

## Scope and files

- `src/oncorna/final_fit_export.py` defines a metadata plan, strict NumPy header inspection, exact-row extraction, a versioned manifest, and sealed-attempt verification. Its only runnable exporter is `export_synthetic`: it requires `SYN-` patient/sample identities and refuses source or output paths under this repository's `data/` directory. There is deliberately no real-export CLI or authorization switch.
- `src/oncorna/final_test.py` now accepts only the reviewed `final_fit_v1` manifest with its completion state and seal. The real-execution and separate test-access guards are unchanged.
- `tests/test_final_fit_export.py` exercises synthetic input only. `tests/test_final_test.py` now creates its fitting-only fixture through the exporter.
- This report records the contract and remaining authorization boundary. Scientific configuration and frozen split were not edited.

## `final_fit_v1` contract

The manifest is local-only JSON in `data/processed/final_fit_v1/manifest.json`, beside `fit_expression.npy`, `success.json`, and `export_state.jsonl`. An exact field set is mandatory; unknown, missing, or incompatible versions fail. The fields are:

| Group | Required fields |
|---|---|
| Contract and method | `version`, `export_method`, `command`, `software` (`exporter_version`, Python, NumPy) |
| Source provenance | `source_expression_sha256`, `source_all_gene_expression_sha256`, `source_schema_sha256`, `frozen_split_sha256`, `approved_cohort` (`identifier`, `path`, `sha256`) |
| Identity and schema | `ordered_fit_identity_sha256`, `ordered_gene_list_sha256`, `fit_patient_count`, `gene_count`, `shape`, `dtype`, `dtype_descr`, `byte_order`, `orientation`, `expression_units` |
| Source layout | `source_array_shape`, `source_array_data_offset`, `source_array_layout` |
| Extraction and output | `selected_payload_sha256`, `repeat_selected_payload_sha256`, `fit_expression_sha256`, `source_hash_verified_at_export`, `source_hash_verification_stage` |

The approved real population is 675 train-plus-validation patients, 169 excluded test patients, and 20,530 original genes. Fitting rows follow the frozen train mapping then validation mapping; columns follow the original all-gene schema order. The artifact is C-order, little-endian float64, shape `[675, 20530]`, orientation `samples_by_genes`, with original `log2(normalized_count + 1)` values. Ordered identities use the existing `identity_digest` encoding; ordered genes use compact UTF-8 JSON. The source-expression and source-array hashes are the frozen recorded values. **The full source-array hash is not recomputed at export**, because doing so would read test-containing bytes; the manifest must state `source_hash_verified_at_export: false` and `source_hash_verification_stage: authorized_test_access`.

## Exact-range algorithm and safety

The exporter verifies the frozen split and schema metadata hashes, partition counts and mappings, unique patient/sample/gene IDs, and the authorized sample set before opening expression data. The NumPy reader opens the source unbuffered, validates the magic, v1/v2 header, exact keys, declared shape, `<f8` dtype, C order, 16-byte data alignment, and exact file length. It rejects unsupported layouts and malformed files. Each fitting sample ID is mapped to one unique in-bounds source row; the fitting set must be disjoint from the test set. It then seeks to `data_offset + row_index * 20530 * 8` and requests exactly one fitting row. It writes those original bytes directly into a new NumPy artifact, without numeric conversion, scaling, filtering, imputation, or model fitting. A second fitting-only read computes an independent selection digest. The output payload digest is checked again during sealed-artifact verification.

These are **application-level byte-range guarantees**: the Python code requests and decodes only authorized fitting rows. The operating system, filesystem cache, storage controller, or hardware may physically read adjacent blocks through read-ahead. This method does not claim physical isolation at those layers. Concurrent source mutation is also outside the current guarantee; the later authorized source checksum and exact-value comparison remain necessary.

The output directory is created exclusively. A pending array is fsynced before promotion; the manifest and success marker are fsynced and atomically placed; the terminal `completed` event is appended only after both exist. Verification requires matching terminal state, manifest and success hashes, the artifact hash, exact manifest fields, valid output layout, and payload digest. An interruption leaves the directory as a visible attempt and a retry with the same path fails. The ignored `data/processed/` location keeps future real artifacts local.

## Synthetic verification

Python 3.11.9: **104 tests passed**, with two existing scikit-learn deprecation warnings. Synthetic tests cover ordered and boundary row selection, absence of test-row requests, invalid layout/header/length, duplicate or wrong membership, every missing manifest field, unknown version/field, altered source and artifact values, byte-deterministic repeat extraction, interruption before and during sealing, and refused retry. `ruff check .` passed. The four touched Python files passed `ruff format --check`. Repository-wide formatting still flags the pre-existing, unrelated `tests/test_download_data.py:35`; it was not edited.

The scientific configuration byte SHA-256 remains `1cadc4bce4d891f904a06409ba12794afa854abf8733276dc4484e2f074baa36`. The frozen split remains unchanged. Real final-test execution remains blocked at the executor, the separate test-access gate, and the real CLI.

## Remaining blockers

1. Review this infrastructure and authorize a separate, controlled real fitting-only export. This checkpoint exposes no real-export entry point.
2. Before any real model fit, verify the produced artifact, its manifest, exact 675/169 membership, and its local-only status. The full original array hash and exact equality to selected original rows remain deferred to the separately authorized test-access stage.
3. Separately authorize any eventual real final-test execution; this checkpoint provides no such authorization and changes no scientific settings.

**Decision:** Ready for human review of the synthetic infrastructure. **No-go** for real artifact creation or final-test execution under this checkpoint.
