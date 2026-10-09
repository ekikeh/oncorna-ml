# Phase 2D-B — final-test infrastructure implementation

**Status:** infrastructure implemented and synthetic-tested on 2026-10-10. The scientific method is frozen for implementation. A real 675-patient fit and access to the 169 test patients' expression or labels remain disabled pending a separate checkpoint. No real final-test model or result was produced.

## Files and frozen configuration

| File | Purpose |
|---|---|
| `.gitattributes` | Keeps the frozen final-test config in LF form on every checkout |
| `configs/final_test_v1.yaml` | Frozen population, model, metrics, reporting, input hashes, and local-only output policy |
| `reports/phase_2d_final_test_protocol.md` | Approved scientific protocol and deterministic reporting rules |
| `src/oncorna/final_test.py` | One-shot state machine, synthetic source, integrity checks, lazy real source, and real-execution gate |
| `scripts/run_final_test_v1.py` | Real-run entry point that currently refuses execution |
| `tests/test_final_test.py` | Synthetic fixtures and safeguard tests only |

The exact byte SHA-256 of `configs/final_test_v1.yaml` is **`1cadc4bce4d891f904a06409ba12794afa854abf8733276dc4484e2f074baa36`**. `.gitattributes` fixes that file's checkout line endings to LF, preserving the hash on Windows. Its real protocol fixes 506 original training plus 169 validation patients for a 675-patient fit, 169 untouched test patients, 20,530 ordered input genes, the strict `>1` expression rule in at least 135 fitting patients, fit-only `StandardScaler`, and L2 logistic regression with `C=1`, `lbfgs`, intercept, no class weighting, `tol=0.0001`, and `max_iter=1000`. Five-class order and macro F1 remain fixed. No scientific parameter was retuned.

## Runner and access boundaries

The one-shot executor validates the protocol and source, performs preflight membership checks, exclusively creates the output directory, and persists `reserved → fitting → fitted → test_access_started → completed` in a flushed state log. A failed stage appends `failed` and preserves a failure record and any partial artifacts. Preflight failures reserve a failure-only directory after the failed check so the attempt remains visible. An existing directory is never overwritten or deleted; a second attempt fails with `FileExistsError`.

The fit stage checks exact fitting patient/sample pairs, ordered gene schema, finite fitting values, filter sample IDs and prevalence count, scaler fit count, and classifier classes. It records warnings and stops on `ConvergenceWarning` before test access. The separate test stage checks the frozen test membership and schema, predicts one batch, aligns probabilities using the estimator's actual `classes_`, validates row sums and class labels, writes the patient-level prediction table, then re-reads and hashes it. Only after that seal does it request labels and calculate fixed metrics. Raw confusion rows represent true classes and columns predictions. Row and column totals must reconcile with support, predicted counts, and the test total.

`SyntheticFixtureSource` accepts only `SYN-` patient and sample IDs. The public executor accepts only that exact synthetic source type in Phase 2D-B. The CLI has no synthetic switch; it requires a separate local authorization record bound to the clean source revision, configuration hash, frozen split hash, and output directory. **Even a matching record cannot start a real fit in this revision:** an additional explicit code-level guard raises `AuthorizationError` before constructing or opening a real source. A second guard remains at the test-access boundary. These guards require a reviewed later code revision and explicit human authorization before any real execution.

The configured final output directory is `data/processed/final_test_v1`, which `.gitignore` excludes. Existing Phase 2B and 2C output directories are never write targets. The frozen split is only a read-only input to the future real source; no split generation path exists in this runner.

## Provenance and reproducibility

The synthetic execution manifest records source revision, configuration hash, input hashes, fitting identity digest, ordered retained-gene digest, model checksum, Python and library versions, solver iterations, warnings, fit and elapsed times, and output checksums. The real source is designed to verify the configured hashes of the source expression, cohort, all-gene array and schema, preprocessing manifest, frozen split, Phase 2B/2C manifests, and their configuration files before loading fitting rows. That real path was deliberately not exercised in this checkpoint.

Historical `configs/default.yaml` provenance uses the existing `verify_split_provenance` path. Its `verify_historical_config` helper accepts exact raw bytes or a provable LF/CRLF-only variant, records the compatibility mode and both byte hashes, and does not replace the historical recorded digest. The new runner carries that compatibility mode into its future manifest. The final-test config itself is bound to its exact byte hash.

## Synthetic verification

Eight focused tests cover disjoint and exact membership, scientific-setting locks, the 20% prevalence boundary and strict `>1` rule, fit-only filter/scaler checks, ordered probability columns and normalization, fixed confusion orientation and `zero_division=0`, output exclusivity and rerun refusal, preflight and fitting failures before test access, failure after test access, persistent state logs, rejection of real execution before source access, synthetic ID enforcement, and identical predictions and metrics from identical synthetic inputs. All fixture identifiers and expression arrays were generated in memory with `SYN-` IDs.

On Python **3.11.9**, the complete suite passed: **85 passed, 2 existing scikit-learn deprecation warnings**. `ruff check .` passed. `ruff format --check .` found one pre-existing formatting difference in `tests/test_download_data.py` line 35; it is outside this phase and was left untouched. The new runner, CLI, and synthetic test file pass Ruff formatting.

## Review boundary and unresolved items

No real fit, test-expression read, test-label read, final prediction, or final metric was run. The real source and real authorization-record path remain unexecuted, so their operational behavior must be reviewed at the later authorized checkpoint. The code-level real-run block is intentional and must not be removed under this Phase 2D-B approval. No patient-level output, model binary, source dataset, credential, or frozen split artifact belongs in this commit.

Human review is requested for the committed implementation and exact config hash. A later checkpoint must decide separately whether to authorize a single real final-test attempt.
