# Phase 2D-F29 — Production fitting activation infrastructure

This checkpoint enables a distinct `production` fitting scope in code. It
does not create an authorization record, claim, model artifact or approval.
The only executable production path requires a later human-approved,
revision-bound local authorization record and an exclusive attempt claim.

## Trust boundary

`run_real_fit` accepts only the reviewed checkout. Preparation independently
derives the 675 fitting patient/sample identities, 169 excluded test patients,
20,530 genes and true PAM50 labels from hash-verified frozen fitting-side
metadata. It checks the exact frozen configuration, split and fitting-artifact
checksums. An authorization record must match the clean executable Git
revision, fixed canonical checkout/output paths, ordered identities and genes,
ordered true-label digest, label-table hashes, frozen model settings, UUIDv4
attempt ID and identifiable human approval reference.

The local claim is created exclusively before fitting-value loading or output
creation. The shared fitting engine independently rechecks the authorization
record, claim bytes, existing Windows DACL and file identity, executable
revision, frozen binding and fitting values before output creation, before and
after fitting, and during sealing. Its first check records authorization,
claim and expression-file identity and modification metadata; later checks
reject changes. Caller-provided callbacks cannot replace these checks.

`synthetic`, `claimed_synthetic`, `production_fixture` and `production` have
distinct manifest scopes. Synthetic scopes cannot be relabeled into production.
The independent verifier accepts `production` only if its authorization,
claim and frozen real fitting binding reconcile. The manifest contains the
authorization checksum and attempt ID; the success marker records the manifest
checksum, and terminal history records both checksums. These locally editable
records establish internal consistency, not independent authentication of a
coordinated rewrite. Verification hashes model bytes without deserializing
joblib. Failed or uncertain outputs and claims remain for manual
reconciliation; the output is never overwritten or retried automatically.

## Scientific and privacy boundary

The fitter retains the frozen prevalence rule (expression strictly greater
than 1 in at least 135 of 675 fitting samples), fitting-only StandardScaler,
and L2 logistic regression with C=1, lbfgs, intercept, no class weighting,
tol=0.0001 and max_iter=1000. No full-cohort or test expression/label source is
used by fitting preparation. The final-test executor and test-access guards
are unchanged. The original full-source checksum and exact fitting-row
comparison remain deferred to separately authorized test access.

## Operating limits

The local approval reference is identifiable but not cryptographically
authenticated. The Windows claim supplies exclusive creation and a restrictive
ACL; it does not prove directory-entry persistence after power loss or provide
a globally enforced one-attempt service. A crash or uncertain state consumes
the attempt until manual review and a new explicit authorization. Python code
running with the operator's permissions is not a security sandbox. A later
test stage that needs to authenticate the saved model against coordinated
rewriting must compare its digest with evidence retained independently of the
editable fitting-output bundle before loading model bytes.

All execution tests for this checkpoint use synthetic identities and values.
The `production_fixture` scope exercises claim-to-seal behavior without
producing a production seal. The real 675-patient model and held-out test set
were not run. Publication and any real fit require separate human decisions.
