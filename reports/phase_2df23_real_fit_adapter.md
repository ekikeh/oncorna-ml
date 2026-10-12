# Phase 2D-F23 — Real fitting adapter integration

## Scope and boundary

This checkpoint adds metadata-only preparation for the final 675-patient fit. It does not execute the real fitting pipeline. `run_real_fit` still raises unconditionally. The final-test executor, CLI, and test-access guards were not changed. No real fitting authorization or attempt claim was created.

## Adapter

`final_model_fit_real.py` derives the ordered fitting identities from the hash-verified frozen split and the ordered genes from the hash-verified all-gene schema. `RealSource.preflight` hashes and verifies the sealed fitting-only expression artifact and fitting-side Phase 2B/2C tables; it skips the full-cohort expression and cohort files. Thus a future real preparation would read fitting-only expression bytes for checksum verification, but no test expression. `RealSource._fit_labels` reads only `true_label` for frozen train and validation pairs. The adapter requires 506 train, 169 validation, 169 excluded test, 20,530 genes, and class counts of 337, 153, 113, 54, and 18 in the frozen class order. It never accepts caller-provided identities, genes, or labels in its public real binding function.

The existing fitting pipeline and sealed output contract remain intact. Its prevalence rule is expression strictly greater than 1 in at least 135 of 675 samples, followed by fitting-only `StandardScaler` and the frozen C=1 L2 logistic regression. The independent verifier hashes sealed files and checks the success marker and state history without loading untrusted joblib.

## Authorization and attempt claim

The old real authorization verifier accepted a caller-supplied expected record. It now derives the expected record from verified metadata, compares the exact `final_model_fit_local_audit_v1` schema, and checks the clean Git revision, canonical checkout and output, frozen hashes, ordered identity/gene/true-label digests, and both label-table hashes. The attempt ID must be UUIDv4 and a human approval reference is required.

The preparation module has a separate claim function. It rechecks the authorization checksum and uses the reviewed Windows `CreateFileW(CREATE_NEW)` claim helper with a creation-time restrictive DACL, reread, and identity verification. A preexisting or interrupted claim blocks another attempt. It creates no model output and does not call the real fitting pipeline. It does not establish power-loss durability or a globally enforced one-attempt rule. A future activation must reconcile the claim and revalidate input stability before accessing fitting values.

## Verification

Synthetic tests cover frozen membership and label derivation, predicted-label exclusion, changed label-table/source metadata, changed expression artifact, authorization mismatch, redirected path, dirty checkout, exclusive claim, simulated ACL failure, interruption and no retry, and locked real fitting/final-test guards. The synthetic fitting tests also cover fitted preprocessing, model sealing, and tamper rejection.

Frozen scientific configuration SHA-256: `1cadc4bce4d891f904a06409ba12794afa854abf8733276dc4484e2f074baa36`.

Frozen split SHA-256: `33e464aa58a69bed96740ef76b011ad23689acf2ed3baef79729b62057a2e26d`.

The original full-source checksum and exact source-row equality remain deferred to separately authorized test access. This checkpoint makes no claim that either was freshly verified. Human approval and a further reviewed activation change are required before any real model fit; held-out test access requires separate authorization after fitted-model sealing.
