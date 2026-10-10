# Phase 2D-F4 — harden the real-export boundary

**Scope:** implementation and synthetic verification only. No real fitting artifact, 675-patient fit, test expression/label access, or final-test run occurred. The real-export function and final-test execution guards remain disabled.

## Changes

- The exact-range exporter now opens the NumPy source once for header validation and both selected-row passes. It records `st_dev`, `st_ino`, size, nanosecond modification time, and nanosecond change/creation time where the platform exposes them. It compares descriptor and path metadata before and after the header and each extraction pass. A mismatch fails the exclusive attempt, retains its state/failure records, and cannot be retried at the same path. The double selected-row digest remains. Frozen split and schema JSON are parsed from the same bytes whose SHA-256 was checked.
- The `final_fit_v1` manifest field set is unchanged. Validation now rejects wrong JSON types recursively, including Boolean-as-integer and integral floats, and requires lowercase 64-character SHA-256 fields. The original source data offset must be an aligned integer. Existing sealed-artifact checks remain.
- The final-test adapter now checks that the fitting manifest has not changed since preflight and independently inspects the original array header and recorded data offset **only after authorized full-source provenance verification begins**. The real final-test executor, separate test-access gate, and real CLI remain locked.
- `final_fit_real.py` defines a metadata-only preparation path. It fixes the config and split hashes; source/schema paths; 675 fitting, 169 excluded test patients, and 20,530 genes; frozen mapping and gene order; and the sole local output path `data/processed/final_fit_v1`. It requires a clean exact Git revision and an exact local authorization record at `data/processed/final_fit_v1_authorization.json`. The record binds the revision, hashes, ordered identity and gene digests, source/schema identifiers, output path, and manifest version. `run_real_export` unconditionally refuses execution. No real-export CLI is enabled and no authorization record was created.

## Synthetic evidence

The full Python 3.11 suite passed with **111 tests** and two existing scikit-learn deprecation warnings. New cases exercise timestamp changes between passes, simulated path replacement and size changes, wrong manifest types, missing or mismatched provenance, wrong fitting counts, incorrect source data offset, and denied real-export invocation. Prior tests continue to cover malformed headers, truncated payloads, authorized row requests, interruptions, retries, and final-test guards. Ruff lint and formatting pass on changed Python files. Repository-wide formatting still flags the unrelated pre-existing `tests/test_download_data.py:35` issue.

## Remaining limits and next authorization

The current function cannot create a real artifact. A later reviewed revision would need to connect the authorization gate to the exact-range writer and receive separate human approval for that exact revision and one fixed output path. An existing output directory, including a failed attempt, must remain a hard stop.

File identity and timestamp checks reduce source-swap and mutation risk; they cannot prove the full original array checksum without reading test-containing bytes. On some filesystems, identity and timestamp precision is limited, and the OS or storage hardware may read ahead of requested fitting byte ranges. The manifest continues to state that the full source hash was **not** verified at export. Full-source cryptographic verification and exact fitting-row comparison remain deferred to the separately authorized test-access stage. No scientific settings or frozen memberships were changed.

**Decision:** ready for independent human review of the still-locked implementation; **no-go** for real artifact creation, real model fitting, or final-test execution at this checkpoint.
