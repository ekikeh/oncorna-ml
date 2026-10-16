# Phase 3B — test-only evaluation implementation

**Status:** implementation and synthetic verification only. This checkpoint does not
authorize access to the 169 held-out patients or execution of the real test.

## Approved fit binding

The test consumer fixes the operator-confirmed SHA-256 values in reviewed source:

| Artifact | SHA-256 |
| --- | --- |
| `data/processed/final_model_fit_v1/pipeline.joblib` | `ea54f2254e72cbd865c737389a500406cb4ced52621c3f64c90176875a967f90` |
| `data/processed/final_model_fit_v1/manifest.json` | `151a40983477cdf63348fef4b70c1b2b874df9c49b5541017fa3087afbb88d53` |

The fitting revision is `aa99b56653ccce0f83ae33bb5af2f454edfe9af5`,
attempt `6dd22770-3986-4d9e-a1f7-bc67f7041c26`, with fitting authorization
SHA-256 `73b8a12c2fa8f1ca2c7ac4e2c331bad6bb9d555c99e499497a9ea9fc8345753e`.
`verify_approved_fit` first hashes both files, then verifies the fitting
authorization, protected attempt claim, sealed outputs, 17,623-gene list,
success marker and terminal state. Historical revision verification permits the
evaluation runner to reside at a new reviewed commit without weakening the
fit's artifact or claim checks. Deserialization happens only after the new
test claim and original-source provenance gate.

## Future authorization contract

The operator must separately approve one exact evaluation revision and one
UUIDv4 test attempt after reviewing this implementation. The local record at
`data/processed/final_test_v1_authorization.json` must match
`expected_test_authorization` exactly: schema and intent; attempt ID and human
approval reference; executable commit; canonical checkout and output path;
frozen config and split hashes; and all five fixed fit identities above. The
record is a future approval artifact, not produced by this implementation.
The runner refuses a dirty checkout, redirected record, existing claim or
existing output directory. It creates and verifies the exclusive protected
claim at `data/processed/final_test_v1_attempt.json` before any test access.
An incomplete claim or output remains for manual reconciliation and blocks
reuse. Neither file has been created by this checkpoint.

## One-time execution order

1. Recalculate trusted model and manifest hashes; verify the sealed fit and
   fitting-only metadata; validate the separate test approval and claim.
2. Create `data/processed/final_test_v1/` exclusively and record `reserved`.
   Preflight reads only metadata and fitting-only artifacts.
3. Record `test_access_started`. Authenticate the original source expression
   checksum and full-cohort array checksum; validate array header, shape,
   offset and ordered gene schema; compare every one of the 675 exported fit
   rows with the authenticated original. Any discrepancy fails before model
   loading or test predictions.
4. Rehash and load the single approved pipeline. Verify its saved filter,
   scaler and classifier against the sealed fit, then select the 169 frozen
   test rows in split order with all 20,530 input genes in original order.
   Invoke `predict` and `predict_proba` once each. No fitting method is called.
5. Write ordered predictions and five probabilities aligned from the saved
   classifier classes. Verify the table independently, write its checksum to
   `prediction_seal.json`, and record `predictions_sealed` with both table and
   seal hashes. Recheck the seal before opening the cohort label file.
6. Only after sealing, verify the cohort checksum and historical split
   provenance, then load the 169 labels. Recheck the predictions. Calculate
   only the fixed metrics and verify Normal-like support of five. Write
   metrics and confusion outputs, a manifest of output hashes, a success
   marker and terminal `completed` state. Independently reconcile the seal,
   prediction schema, probability rows, manifest and state history.

The frozen primary metric is five-class macro F1. Supplementary outputs are
balanced accuracy, overall accuracy, class precision/recall/F1/support,
raw and row-normalized confusion matrices, predicted class counts and aligned
uncalibrated probabilities. The Normal-like test support is five, so one case
changes its recall by 20 percentage points. No model, threshold, calibration
or feature adjustment is permitted after inspecting the test.

All patient-level outputs remain local. A failed or uncertain attempt leaves
its failure record, state log and partial files intact. There is no automatic
retry, replacement or second evaluation; any recovery requires a documented
root cause, a separately reviewed version and fresh human authorization.
