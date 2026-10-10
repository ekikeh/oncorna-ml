# Phase 2D-F6 — locked real-export activation path

**Scope:** code integration and synthetic verification only. No real fitting artifact, real 675-patient fit, held-out expression or label access, or final-test evaluation occurred.

## Integration and authorization

`src/oncorna/final_fit_export.py` now has one private exact-range writer used by the synthetic entry point and referenced by the real entry point. The writer's `real=True` branch unconditionally raises **before directory creation or source access**. The synthetic branch continues to require `SYN-` identities and source/output locations outside repository `data/`. A later reviewed code revision must explicitly change this branch; a local authorization record alone cannot activate it.

`src/oncorna/final_fit_real.py` now connects `run_real_export` to `prepare_real_export` and then to the shared writer. Preparation requires the one fixed, Git-ignored authorization-record path, exact current Git revision and clean tree, frozen configuration and split checksums, recorded array and schema identifiers/hashes, ordered fitting identity and gene digests, `final_fit_v1`, and the sole output path `data/processed/final_fit_v1`. It derives 506 train plus 169 validation identities, excludes 169 test identities, and retains all 20,530 original genes from verified frozen metadata. No caller-supplied membership list enters this path. A pre-existing output directory is rejected; the writer creates new outputs exclusively and preserves failed attempts.

The manifest field set is unchanged. Its source-array and raw-expression hashes remain recorded provenance, **not freshly computed full-source hashes**. The extraction computes fitting-row and output checksums, repeats fitting-row reads, and verifies the sealed manifest, success marker, and artifact. The source descriptor and path metadata are checked during extraction. Full-source checksum and exact source-row equality remain deferred to a separately authorized test-access stage. The final-test executor, real CLI, and separate test-access gate were not changed.

## Synthetic evidence and remaining boundary

Tests build an 844-identity, 20,530-gene **synthetic metadata-only** fixture. Preparation succeeds without an all-gene expression file; real-mode execution then refuses before creating output. Synthetic tests reject wrong revision, dirty tree, wrong record fields and output path, bad membership, and direct real-mode invocation. Existing synthetic tests cover exact authorized row requests, untouched test rows, source mutation/replacement, unsupported NumPy layouts, interruptions, retry refusal, and sealed output verification.

The full Python 3.11 suite passed with **115 tests** and two existing scikit-learn deprecation warnings. Ruff lint and formatting pass on changed Python files. Repository-wide formatting still flags the unrelated pre-existing `tests/test_download_data.py:35` issue. The frozen configuration and split SHA-256 values remain `1cadc4bce4d891f904a06409ba12794afa854abf8733276dc4484e2f074baa36` and `33e464aa58a69bed96740ef76b011ad23689acf2ed3baef79729b62057a2e26d`.

A later checkpoint must review a new exact Git revision that changes the unconditional code guard, verify that the fixed authorization path is the only real entry, and obtain explicit human authorization for one attempt at the fixed output path. File metadata checks do not replace full cryptographic source verification, and OS/storage read-ahead may physically touch adjacent blocks. **Real export, model fitting, and final-test evaluation remain unauthorized and disabled in this revision.**
