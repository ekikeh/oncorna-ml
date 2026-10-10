# Phase 2D-F7 — F6 publication and local activation review

**F6 publication:** The reviewed locked F6 commit `1cba25528bd445d2e475142f2b4616955820c664` was pushed by a non-force fast-forward from `42cafb9695a3d6834e1b56540eb04f143c9c7930`. An independent `git ls-remote` read returned the exact F6 SHA for GitHub `main`. Its working tree was clean; tests, lint, four-file inventory, Git ignore rules, and real-export/final-test guards were verified before push.

## Local activation patch

The shared exact-range writer now accepts real mode only after it independently re-runs `prepare_real_export` against the fixed authorization path and requires the resulting plan to equal the caller's plan. This check rebinds the reviewed Git revision, clean tree, frozen configuration and split hashes, recorded source/schema identifiers and hashes, ordered 675 fitting identities, 169 excluded test identities, 20,530 genes in schema order, manifest version, and fixed local output path. The real entry point prepares the plan and passes it to that writer; callers without the same validated plan and authorization cannot reach output creation.

Immediately before output creation, real mode rechecks the current revision, clean tree, fixed authorization path, and exact authorization-record fields. It exclusively creates `data/processed/final_fit_v1_attempt.json` as a persistent claim containing the hash of those checked record bytes, ordered identity and gene digests, output path, and claim time. The fixed marker and authorization record are Git-ignored. An existing marker or output directory blocks another attempt; an interruption after claiming leaves the marker visible even if output creation fails. The writer then uses the already-reviewed header, source-stability, fitting-row-only, byte-copy, exclusive output, state-log, manifest and success-seal logic. Its `final_fit_v1` schema and source-hash deferral are unchanged. The full original source checksum is recorded but **not freshly verified** at export. No filtering, transformation, scaling, imputation, model fit, or test prediction was added.

The activation patch changes only `src/oncorna/final_fit_export.py`, `src/oncorna/final_fit_real.py`, synthetic tests, and this report. It does not touch `src/oncorna/final_test.py`; the real final-test executor, test-access gate, and real CLI remain locked. A fitting-export authorization record does not match the separate final-test authorization schema.

## Synthetic verification

The full Python 3.11 suite passed: **119 tests**, with two existing scikit-learn deprecation warnings. Synthetic tests exercise authorized fitting-row extraction and sealing, absence of test-row requests, revision and dirty-tree rejection, missing/mismatched authorization, wrong membership and output paths, existing output/attempt refusal, source mutation or replacement, malformed NumPy layouts, interrupted attempts, and rejection of the fitting-export record by the final-test authorization checker. Ruff lint and changed-file formatting pass. Repository-wide formatting still flags only the pre-existing unrelated `tests/test_download_data.py:35` issue.

The frozen configuration SHA-256 remains `1cadc4bce4d891f904a06409ba12794afa854abf8733276dc4484e2f074baa36`; the frozen split SHA-256 remains `33e464aa58a69bed96740ef76b011ad23689acf2ed3baef79729b62057a2e26d`. No real fitting artifact or local authorization record was created. No real expression values or test labels were accessed.

## Human approval boundary

This patch is for local review only and is not pushed. It makes a real fitting-only export **technically executable only if** a future local authorization record names this exact export-enabled commit and its fixed inputs/output. No such record exists. A human must independently review the committed diff and explicitly authorize one export attempt before that record is created or the real entry point is run. Final-model fitting and final-test execution require separate approvals. File identity checks do not replace the deferred full-source checksum, and physical storage read-ahead remains possible.

**Recommendation:** GO to review the local activation commit; **NO-GO to run a real export now**. Stop before any real fitting artifact, model fit, or test evaluation.
