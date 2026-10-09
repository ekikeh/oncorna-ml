# Phase 2D-D — final-test safety corrections

**Status:** implementation and synthetic verification only, 2026-10-10. The real final-test guards remain in place. No real 675-patient fit, real test expression or label access, final prediction, or final evaluation occurred. Nothing was pushed.

## Corrections

| Phase 2D-C finding | Correction |
|---|---|
| Preflight read full-cohort bytes | Preflight now hashes metadata and fitting-only artifacts. The all-gene array, source expression, cohort file, and historical full matrix are verified only after `test_access_started`. These full-file checks are mandatory and fail before test prediction. |
| Fitting labels parsed test rows | Fitting labels now come from the existing Phase 2B `oof_predictions.tsv` (the selected `C=1.0` rows for 506 training patients) and Phase 2C `validation_predictions.tsv` (169 validation patients). Their hashes are verified against their already-frozen manifests before parsing; exact patient/sample membership and duplicate checks apply. Fitting no longer opens the full cohort label file. |
| `completed` preceded sealing | All result hashes are recorded and verified; the manifest and success marker are written through fsynced temporary files and atomic replacement. The terminal `completed` state is appended only after those checks. `verify_completed_run` requires the terminal state, success marker, manifest digest, and all output hashes to agree. A crash before either seal leaves a visibly incomplete audit state. |

The fitting expression adapter also previously opened the all-patient `.npy` file before test authorization, despite selecting only fitting rows. It now requires `data/processed/final_fit_v1/fit_expression.npy` and an exact provenance manifest containing the original array hash, frozen split hash, ordered fit identity digest, ordered gene digest, shape, dtype, and artifact hash. Preflight and fitting recheck that fitting-only file. After the authorized full-file provenance check, the runner compares the fitting-only values with the corresponding verified original rows before any test prediction. **This real fitting-only artifact does not exist yet and was not created in this checkpoint.**

The artifact must be supplied at a separately reviewed checkpoint from a controlled fitting-only export that does not read test values. Its provenance manifest must be reviewed before any real fit. The later authorized stage independently checks exact value equality against the original array. A manifest assertion alone is not proof of equality before that later check. The frozen configuration and scientific method were not changed.

## Access and audit model

1. Metadata-only preflight checks split counts and disjointness, all-gene schema, metadata hashes, historical configuration LF/CRLF compatibility, existing training/validation prediction hashes, and the fitting-only expression manifest and checksum. It does not open test-containing expression or label files.
2. Fitting reads only the fitting-only expression artifact and the verified training/validation true-label columns. The filter, scaler, and classifier receive exactly the frozen fitting IDs.
3. The executor records `test_access_started` before enabling the source's test-access method. The source then verifies all deferred full-cohort hashes and historical matrix provenance. A mismatch fails the attempt before test prediction. Test expression is unavailable until this succeeds; test labels additionally require saved and re-read predictions.
4. The final manifest hashes the model, retained gene list, training manifest, prediction table, metrics, and confusion matrices. The success marker names the manifest hash; the terminal audit event names both hashes. The run-state log is excluded from the manifest's output hashes to avoid a circular digest. A run is successful only if `verify_completed_run` accepts the terminal state and every seal.

Exclusive output-directory creation and `x` mode for individual artifacts continue to prevent silent overwrite. Failed and interrupted directories remain visible; no automatic retry is provided. The real CLI still refuses execution even with an authorization record, and the executor still rejects a real source before preflight.

## Synthetic evidence and limitations

Synthetic adapter fixtures omit the full-cohort array, source matrix, and cohort label file entirely. Preflight and fitting-only loading succeed from metadata, a fitting-only array, and fitting-only label artifacts. Attempts to verify full provenance before the test-access stage, to load test expression before full verification, or to load test labels before sealed predictions fail. After the synthetic stage begins, missing full-cohort files cause provenance verification to fail rather than being skipped. A changed fitting-only array fails before loading. A separate synthetic test changes held-out expression and labels and obtains the same fitted model checksum and retained-gene digest. Interruptions before manifest creation and after manifest creation but before the success marker both leave a terminal state other than `completed`; the completion verifier rejects both and reruns cannot reuse their directories.

The real fitting-only expression artifact is the remaining preparation blocker. The dormant real adapter and full provenance/equality path were tested only with synthetic inputs or reviewed statically. Existing Phase 2B/C training and validation prediction files were inspected for their label columns and coverage; no real expression values were loaded and no real model was fitted. A later checkpoint must review the fit-only artifact's creation and provenance, then separately authorize any real execution. The Phase 2D-C safety review should be repeated on the resulting commit before publication or final-test approval.

## Validation and Git scope

The frozen `configs/final_test_v1.yaml` byte SHA-256 remains **`1cadc4bce4d891f904a06409ba12794afa854abf8733276dc4484e2f074baa36`**. Phase 2D-D changes only `src/oncorna/final_test.py`, `tests/test_final_test.py`, the Phase 2D protocol report, the Phase 2D-B report inventory/safety notes, and this report. No data, split, model, patient-level prediction, or credential artifact is eligible for this commit.

On Python 3.11.9, the focused synthetic suite has 13 passing tests and the complete suite has **90 passing tests** with two existing scikit-learn deprecation warnings. `ruff check .` passes, and the changed Python files pass Ruff formatting. Repository-wide Ruff formatting retains the pre-existing unrelated `tests/test_download_data.py:35` difference.

**Review recommendation:** GO for independent review of the corrected, still-locked infrastructure. Publication requires explicit human approval after that review. Real final-test execution remains NO-GO.
