# Phase 2D-F27 — Direct-call fitting boundary repair

## Root cause

The F25 `_fit_and_seal` helper accepted caller-provided `checkpoint` and
`validate_manifest` callbacks. Direct callers could substitute no-op callbacks,
provide arbitrary values and provenance, and write a self-consistent model seal
without the real-fitting authorization or Windows attempt claim.

## Correction

- The shared engine now produces only explicit `synthetic` or
  `claimed_synthetic` results. It rejects all production execution scopes. The
  public real-fitting entry point also remains unconditionally disabled.
- Direct synthetic fitting checks synthetic identities, the input artifact and
  output location inside the engine. A synthetic result cannot occupy the fixed
  production output directory or claim a production approval.
- Claimed synthetic fitting independently re-derives the fitting metadata and
  checks the authorization checksum, exact authorization record, exclusive
  claim contents, revision, paths and fitting values before output creation
  and at subsequent fitting and sealing checkpoints. These checks run even
  when caller-supplied callbacks do nothing.
- The manifest records the authorization checksum. The success marker seals
  that manifest, and the terminal state seals both. The independent verifier
  distinguishes synthetic and claimed synthetic scopes, reconciles the claim
  with the frozen fitting binding, and rejects production scope while real
  activation remains locked. The existing output filenames and manifest
  version remain `final_model_fit_v1`; the authorization checksum is a new
  mandatory manifest field before any real model artifact exists.

## Synthetic verification

Python 3.11: 183 passed, one Windows directory-symlink test skipped because
the required privilege is unavailable. Ruff lint passed. Changed-file Ruff
formatting passed. Repository-wide formatting still reports the pre-existing
`tests/test_download_data.py` discrepancy. Regression tests cover direct
no-op callbacks, missing or forged claims, changed values, synthetic output
isolation, coordinated metadata re-sealing, interruptions, and model/manifest/
success-marker tampering.

## Scientific and operational boundary

The configuration, frozen split and fitting expression bytes were unchanged.
No real model fit, authorization record, attempt claim or model output was
created. The final-test executor and access guards were not changed. This
repair is ready for independent publication review; a separate reviewed
activation change and explicit human approval are still required before a
real fit. Python code running with the operator's filesystem permissions is
not a security sandbox; the boundary prevents accidental API bypass and
unsupported trusted seals, not arbitrary malicious process modification.
