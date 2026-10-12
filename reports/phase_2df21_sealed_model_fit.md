# Phase 2D-F21 — sealed final-model fitting infrastructure

## Boundary

`oncorna.final_model_fit` implements a synthetic-only fitting and sealing path.
`run_real_fit` raises unconditionally, before reading a real authorization or
fitting value. The existing final-test runner, CLI, and test-access guards are
unchanged. No F21 code enables real fitting or held-out access.

The future real fitting entry point must derive its fitting identities, gene
order, true labels, and input hashes from the frozen split, schema, sealed
`final_fit_v1` artifact, and hashed Phase 2B/2C `true_label` tables. It must not
accept a caller-selected population. The current synthetic path accepts only
`SYN-` identifiers and paths outside the project's `data` directory.

## Scientific method

The fitting stage calls the existing `make_pipeline` and `modeling_parameters`.
The frozen pipeline applies `GenePrevalenceFilter` (`value > 1` in at least
`ceil(0.20 × n_fit)` samples), then `StandardScaler`, then L2 logistic regression
with `C=1`, `lbfgs`, intercept, no class weight, tolerance `0.0001`, and at most
1,000 iterations. It records the actual retained-gene list and count. It does
not compare candidates or use test or validation scores to alter settings.

## Proposed real authorization

The distinct `final_model_fit_local_audit_v1` schema binds a UUIDv4 attempt,
human approval reference, exact reviewed Git revision, canonical checkout,
fixed `data/processed/final_model_fit_v1` output, frozen config and split
hashes, fitting artifact hash, ordered fitting identities, original gene order,
ordered ground-truth label digest, hashes of both fitting-label tables, frozen
class order and parameters, and manifest version. The authorization is specific
to fitting and cannot authorize test access. `verify_real_fit_authorization`
checks the clean revision, fixed paths, absence of prior output or attempt,
frozen hashes, and exact record. This verification alone cannot start a fit.
Before enabling real execution, a separately reviewed integration must derive
the expected record exclusively from frozen metadata and approved fitting-side
files, create the local attempt claim before output or source access, and bind
that claim to the authorization checksum. F21 does not implement or create a
real claim.

## Output contract and recovery

A synthetic attempt exclusively creates one directory containing
`pipeline.joblib`, `training_gene_list.tsv`, `training_manifest.json`,
`manifest.json`, `success.json`, and `fit_state.jsonl`. The fitting manifest
records input identities and hashes, scientific parameters, software versions,
model and output hashes, retained genes, attempt and authorization references,
and the fact that full-source provenance and exact source-row equality remain
unverified. The training manifest records class counts, filter and scaler fit
counts, retained gene IDs, classifier iterations, warnings, and digests of
fitted numeric state. Completion is recorded only after the model, manifest,
and success seal verify together. Failure evidence and partial files remain;
an existing directory prevents automatic retry.

`verify_sealed_fit` hashes model bytes and the other outputs, validates the
manifest, gene list, success marker, and terminal history without loading the
serialized model. `joblib` files can execute code when deserialized. A later
test stage may load only this locally produced, independently hash-verified
artifact after separate review and authorization; arbitrary external model
files are unsupported.

File flushes and exclusive creation do not prove Windows directory-entry
durability after power loss. A crash or uncertain seal consumes the local
attempt under the proposed policy and requires human reconciliation before any
new approval; the local record is not a globally enforced one-attempt service.

The future test stage must verify this seal, authenticate the full original
source and fitting-row equality at authorized test access, load the saved model
without refitting, seal test predictions before labels, and perform the single
prespecified evaluation. A provenance mismatch must invalidate the run before
prediction or scoring. Neither this stage nor its authorization exists in F21.
