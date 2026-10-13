# Phase 2D-F25 — Controlled final-model fitting integration

## Execution boundary

The public `run_real_fit` entry point calls the new claimed-fit integration, whose real branch raises `AuthorizationError` before reading authorization records or fitting values. This commit does not authorize or enable a real 675-patient fit. No real authorization, attempt claim, or model-output directory was created. The final-test executor, CLI, configuration, and test-access guard are unchanged.

## Frozen input binding and local audit

The real adapter independently derives the ordered 506 training and 169 validation identities, 169 excluded test identities, original 20,530 genes, and five-class ground-truth targets from the frozen split, schema, sealed fitting-only expression artifact, and hash-verified Phase 2B/C `true_label` tables. The exact `final_model_fit_local_audit_v1` record binds these digests, both label-table hashes, scientific settings, canonical checkout/output, exact Git revision, UUIDv4 attempt ID, and human approval reference. Neither caller-supplied labels nor caller-supplied patient or gene ordering is accepted in the real binding.

The integration first verifies authorization and clean revision, then creates the fixed, exclusive Windows local claim with the reviewed creation-time ACL. It compares the claim and authorization again before loading fitting values. It holds one fitting-artifact descriptor open, checks descriptor/path identity, size, high-resolution timestamps and SHA-256, and rechecks the artifact and frozen metadata before fitting, after fitting, before sealing, and after completion. The source array is the 675-row fitting-only artifact; no full-cohort expression or test labels are opened. Checksum verification of the fitting-only artifact occurs during preparation, before the claim; model fitting and output creation occur only after the claim.

The local claim remains after failure or interruption. A preexisting claim or output blocks a repeat. The local audit does not cryptographically authenticate a human approver, globally enforce one attempt, or prove Windows directory-entry persistence after power loss. Uncertain attempts require manual reconciliation and new approval.

## Fitted method and output

Synthetic integration uses the same `_fit_and_seal` implementation as the future real branch. The scientific method remains the fixed prevalence filter (expression strictly greater than 1 in at least 20% of fitting samples, or 135 of 675), fitting-only `StandardScaler`, and C=1 L2 logistic regression with the frozen solver and other settings. No candidate selection or test access occurs.

The common fitter saves the fitted pipeline, retained-gene list, training manifest, provenance manifest, success marker and state history under the fixed local-only output contract. Before sealing, the real integration compares manifest provenance to the approved claim. The independent verifier checks all files, hashes, state history, and success marker without deserializing untrusted joblib. A future test stage must additionally bind this seal to the approved authorization and verified source provenance before loading the locally produced model.

## Verification and remaining gate

Synthetic tests exercise complete claim-to-seal integration, input/claim mutation, same-byte source replacement, interrupted fitting and sealing, refusal to retry, tampered model and manifest detection, and the unchanged final-test guard. The synthetic fixture has no full-cohort expression or test-label files, demonstrating that neither is needed by this path.

Frozen configuration SHA-256: `1cadc4bce4d891f904a06409ba12794afa854abf8733276dc4484e2f074baa36`.

Frozen split SHA-256: `33e464aa58a69bed96740ef76b011ad23689acf2ed3baef79729b62057a2e26d`.

Fitting-only artifact SHA-256: `54db5e3d8d931c22c19c5ac871a8da720f06c8f4f921759cb2d2031f3ae2382a`.

The real branch remains blocked by an unconditional guard. A separate human decision and reviewed code revision must remove that guard, followed by a new authorization record bound to that exact revision. The original full-source checksum and exact fitting-row comparison remain deferred to separately authorized test access. No final-model fitting or held-out evaluation is authorized by this checkpoint.
