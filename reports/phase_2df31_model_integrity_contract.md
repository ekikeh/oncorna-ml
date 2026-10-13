# Phase 2D-F31 — Local final-model integrity contract

The F30 finding is reproducible with a synthetic `production_fixture`: changing
the saved model bytes alone fails verification, while changing those bytes and
updating the manifest, success marker, and terminal-state hashes together can
produce another internally consistent bundle. Those records are all editable
by the local operator. The self-contained verifier does not authenticate that
bundle against a coordinated rewrite.

`verify_sealed_fit` checks the required files, exact manifest schema, model
and other artifact hashes, fitting metadata and fixed scientific settings,
authorization and claim agreement for claimed scopes, absence of failure or
pending files, and the completed state sequence. It hashes model bytes without
deserializing joblib. Missing files, isolated modifications, contradictory
claims, incomplete attempts, and synthetic-to-production scope promotion fail.
An internally consistent rewrite of all editable output hashes and records can
still pass. The verifier establishes local consistency and provenance-field
agreement, not an independently anchored model identity.

Before a later test stage loads the model, any requirement to authenticate its
identity against coordinated rewriting needs a model or bundle digest retained
independently of `data/processed/final_model_fit_v1/` and compared with the
current bytes. Merely adding another checksum inside that directory would not
provide this guarantee. No signing, ledger, final-test access, or new approval
mechanism is introduced here.

The production fitting authorization, exclusive Windows claim, frozen input
checks, direct-call boundary, and completion checks remain in force. This
clarification changes no scientific setting, split membership, fitting-only
expression artifact, or final-test guard. A real fitting attempt and held-out
test access each still require separate human authorization.
