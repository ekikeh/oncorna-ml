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
contains no secrets. The third-party preview's reuse terms require independent review
for a full repository release; MIT does not grant dataset redistribution rights.

Cite [CITATION.cff](../CITATION.cff) and credit TCGA/UCSC Xena using the
[provenance and attribution](data_provenance.md). **Research and education only;
not for clinical use.**

## Steps after separate release authorization

1. Independently review the candidate commit, scientific scope, preview attribution,
   archive inventories and recorded hashes; resolve any redistribution concerns.
2. If authorized, push the reviewed changes and require both platform CI jobs to pass
   on the exact intended release commit. Do not substitute the baseline CI result.
3. Set the actual release date in CITATION.cff when known, review and commit it,
   rerun affected checks and both CI jobs, then rebuild/audit the final archives.
   Record new hashes because citation changes alter archive contents.
4. Obtain approval for that exact commit and artifacts. Only then create and push
   the authorized `v0.1.0` tag and create the GitHub release with these scope notes
   and only audited attachments plus checksums. Inspect the auto-generated source
   archives as full-repository distributions, including preview licensing.
5. Package-index or container publication needs its own explicit authorization.
   No release operation is authorized by successful local verification.
