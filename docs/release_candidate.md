# v0.1.0 release-candidate notes

Status: proposed first research release; not tagged or published. Preparation date:
2026-10-27. This is not a historical release date. Citation version remains 0.1.0;
`date-released` is omitted until an actual release is authorized and made.

## Contents and attribution

The public repository contains historical research source, methodology, approved aggregate
results, synthetic/mocked tests and reproducibility infrastructure. Its attributed real
training preview is for structural inspection; upstream TCGA/Xena terms apply independently
of the MIT code license. Review [provenance](data_provenance.md) before redistribution.

The installable wheel contains only `oncorna/__init__.py`, `prediction.py`,
`synthetic_demo.py` and standard metadata/license files. The 31-file source archive
contains those modules and the explicitly reviewed build/demo support files. It is not
the full research checkout. The real preview, full patient data, fitted pipeline,
held-out predictions, authorization records, claims, seals and private outputs are excluded.
The release notes and verification report remain repository documentation, outside the
unchanged distribution allowlist. See [release scope](release_scope.md).

## Verified scientific scope

This study reconstructs five existing expression-derived PAM50 labels in TCGA-BRCA.
There are 844 approved samples: 506 training, 169 validation and 169 held-out test.
Final fitting used 675 training-plus-validation samples. The ordered input has 20,530
genes; the validation fit retained 17,613 and the distinct final fit retained 17,623.
Development used four-fold CV, not nested CV. The approved one-time held-out results
are macro F1 **0.7934**, balanced accuracy **0.7556**, accuracy **0.8521**
(**144/169** correct), with Normal-like recall **2/5** among five PAM50 classes.
These numbers are transcribed from the [approved aggregate report](../reports/phase_3d_aggregate_results.md),
not recomputed from protected inputs. The held-out evaluation remains permanently closed.

Outputs are uncalibrated. There is no external validation or demonstrated clinical
validity. Expression-derived labels share signal with predictors; upstream normalization
cannot be isolated across partitions and the prevalence threshold arose in cohort-wide
QC. An entirely leakage-free preprocessing history is not established. No DESeq2,
GSEA, SHAP, nested CV, calibration or survival analysis is claimed as completed.
This is not biomarker discovery, a validated diagnostic, a clinical risk predictor or
a prognostic survival model. See the [model card](model_card.md).

## Install and run the offline synthetic demonstration

From a reviewed checkout, use CPython 3.11.17 on glibc Linux x86_64:

```bash
python3.11 -m venv .venv-demo
.venv-demo/bin/python -m pip --isolated install --require-hashes --no-deps -r requirements/runtime-linux-py311.lock
.venv-demo/bin/python -m pip --isolated install --no-deps dist/oncorna_ml-0.1.0-py3-none-any.whl
.venv-demo/bin/python -m pip check
.venv-demo/bin/python -B -m oncorna.synthetic_demo
```

On Windows AMD64 use CPython 3.11.9, `py -3.11 -m venv .venv-demo`,
`.venv-demo\Scripts\python.exe` and `requirements/runtime-windows-py311.lock`.
The wheel must first be built with the development lock as described in
[packaging](packaging.md); it is not available from a published release yet.
The CLI reports `synthetic_demo_only`: 100 training / 40 evaluation samples,
12 invented genes, five classes, default seed 20261026. Synthetic performance is
not research performance. No real fitted model is supplied or loaded.

Only the two documented 64-bit profiles are validated. macOS, ARM and musl are
unsupported by these locks. Historical explicit-logistic-penalty deprecation warnings
remain; the new locks do not establish historical scientific bitwise reproducibility.
Windows installations may need pip's documented system-certificate trust option;
TLS verification must remain enabled. See [packaging](packaging.md) for details and
the restricted, synthetic-only Linux container commands.

## Evidence and release controls

The clean starting revision and published main were both
`496aeeab552662e7572aa55ecf73d4e9711326c0`. GitHub Actions
[run 40379479194](https://github.com/ekikeh/oncorna-ml/actions/runs/40379479194)
completed successfully at that revision for Ubuntu 24.04 / Python 3.11.17 and
Windows 2025 / Python 3.11.9. This is baseline evidence, not CI on the later
local candidate commit. Workflow permissions are read-only, actions are pinned,
checkout credentials are not persisted, and there are no publication steps.

The [Phase 5A verification report](../reports/phase_5a_release_candidate.md) records
candidate archive hashes, tests and limitations. Exact allowlists and content audits
apply to both archives. Scanners supplement review and cannot prove arbitrary text
contains no secrets. MIT does not grant third-party dataset redistribution rights;
the separate evidence for this preview is recorded below.

Phase 5C records the [verified TCGA/Xena redistribution evidence](data_provenance.md#verified-redistribution-evidence-reviewed-2026-10-27)
and its scope, attribution and non-reidentification conditions. Preview identifiers
are linkable, not anonymous. This supports inclusion of the documented preview in
repository snapshots without extending MIT to datasets or permissions to other resources.
The preparation baseline `6ccfc76156933c7b5564547c762f25994e79df16` passed both
hosted platforms in [run 40390676509](https://github.com/ekikeh/oncorna-ml/actions/runs/40390676509):
Linux 387 passed / 13 skipped / 27 warnings; Windows 336 passed / 64 skipped /
2 warnings. That run does not cover later provenance-preparation changes.

Cite [CITATION.cff](../CITATION.cff) and credit TCGA/UCSC Xena using the
[provenance and attribution](data_provenance.md). **Research and education only;
not for clinical use.**

## Steps after separate release authorization

1. Independently review the candidate commit, scientific scope, preview attribution,
   archive inventories and recorded hashes; resolve any redistribution concerns.
2. If authorized, push the reviewed changes and require both platform CI jobs to pass
   on the exact intended release commit. Do not substitute the baseline CI result.
3. Leave `date-released` absent while the release date is unknown. Once the actual
   publication date is established and citation editing is authorized, set it to that
   calendar date (YYYY-MM-DD), review and commit the change. Do not use a preparation
   date as a substitute. If publication is deferred to another date, correct the date
   and repeat the affected gates before release. Record the exact final commit SHA
   and require both hosted CI jobs to pass on that SHA after any citation change.
4. From that clean final commit, build/audit the wheel and source archive with the
   locked tools, verify the clean install and deterministic demo, and compute new
   SHA-256 hashes for the exact files to attach. Preserve the actual audited bytes;
   a later rebuild is not automatically equivalent. Review a full-repository snapshot
   before tagging: it includes the attributed real preview, historical source and tests,
   unlike the 8-file wheel and 31-file package source archive. Check preview provenance
   and absence of protected artifacts. Obtain explicit authorization naming the final
   commit, `v0.1.0` tag and attachments before any release operation. Only then create
   and push that tag at the approved commit and create the GitHub release. Attach only
   the audited wheel, package source archive and their matching checksum manifest;
   verify uploaded hashes and the automatically generated repository snapshots.
5. Package-index or container publication needs its own explicit authorization.
   No release operation is authorized by successful local verification.
