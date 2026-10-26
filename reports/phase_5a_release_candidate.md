# Phase 5A — v0.1.0 local release-candidate verification

Preparation date: 2026-10-27. No release date is asserted.

## Baseline and boundaries

The starting checkout was clean at `496aeeab552662e7572aa55ecf73d4e9711326c0`.
Read-only GitHub API checks confirmed main at the same revision and
[Actions run 40379479194](https://github.com/ekikeh/oncorna-ml/actions/runs/40379479194)
completed successfully on that SHA. Both jobs succeeded:
`checks (ubuntu-24.04, 3.11.17, linux)` and
`checks (windows-2025, 3.11.9, windows)`.
That evidence covers the baseline, not the later local documentation commit.
The workflow retains read-only permissions, immutable action references, disabled
credential persistence, platform hash locks and no publication operations.

All research source, tests, locks, configuration, scientific reports, prediction
contract and synthetic CLI remain unchanged. The three installable module sources
were compared directly with the expected baseline. Review used approved aggregate
reports only; no patient-level inputs, fitted binaries or sealed predictions were
opened. No evaluation, refitting or new research analysis occurred.

## Scientific review

README, model card, release scope, reproducibility appendix and approved Phase 3D
report agree on expression-derived TCGA-BRCA PAM50 reconstruction: 844 samples,
506/169/169 split, final fit on 675, ordered 20,530-gene input, 17,613 retained
validation-fit versus 17,623 final-fit genes, and four-fold development CV rather
than nested CV. Held-out macro F1 is 0.7934, balanced accuracy 0.7556, accuracy
0.8521 (144/169), with Normal-like recall 2/5 across five classes.
These are existing aggregate figures, not newly calculated scientific results.

The documents disclose uncalibrated outputs, no external or clinical validation,
small-class uncertainty, upstream normalization limitations and historical
cohort-wide prevalence-rule selection. They do not claim biomarker discovery,
diagnostic validity, clinical risk prediction or survival prognosis. DESeq2, GSEA,
SHAP, nested CV, calibration and survival analysis remain uncompleted proposals.
The locked synthetic environment does not establish historical bitwise reproducibility.

## Identity and exact changed-file inventory

The proposed version is **v0.1.0**, already consistent in package metadata,
`oncorna.__version__`, citation, archive auditor and CI archive paths. No version
or archive-path change was needed. The premature citation `date-released` was
removed; it must be set to the actual date only for an authorized release.

Exactly four tracked files change:

1. `CITATION.cff`: omit the premature release date.
2. `README.md`: link candidate notes and clarify the installable source subset.
3. `docs/release_candidate.md`: concise release notes, installation, limitations,
   supported profiles, attribution and subsequent release steps.
4. `reports/phase_5a_release_candidate.md`: this verification evidence.

The latter two documents remain outside the installable archive policy. There is
no source, metadata-version, lock, workflow or distribution-policy change.

## Distribution and privacy review

The exact allowlists remain eight wheel files and 31 source-archive files.
The wheel has three demo modules plus standard metadata/license files. The source
archive adds only reviewed demo/build support text, not historical evaluation modules.
The real training preview, full patient data, fitted pipeline, predictions,
authorization records, claims, seals and private outputs are excluded.
MIT covers project code, not third-party datasets. Full-repository release archives
include the existing real preview, whose upstream redistribution terms require
independent review; the demo archives exclude it. Pattern scanners supplement review
and cannot establish absence of every possible secret.

## Candidate validation

Build profile: native Windows AMD64 CPython 3.11.9. A new dedicated development
environment installed all 33 packages from `dev-windows-py311.lock` with
`--require-hashes --no-deps`; editable installation used `--no-build-isolation`
and `pip check` passed. `python -m build --no-isolation` built the source archive
and then the wheel from that archive using the locked build backend.

| Check | Result |
|---|---|
| Full locked Windows synthetic/mocked suite | 335 passed, 65 skipped, 2 known warnings |
| Ruff | Passed |
| Formatting (`src scripts tests tools`) | All 58 files passed |
| Extracted source archive focused tests | 30 passed |
| Exact wheel/source inventory and sensitive-content audits | Passed: 8 / 31 files |
| Safe extraction after audit | Passed for both archives; destinations constrained to dedicated extraction roots |
| Package metadata and citation consistency | Passed: oncorna-ml 0.1.0, Python >=3.11,<3.12, MIT, correct repository |
| TOML / workflow YAML / citation field checks | Passed |
| Local Markdown link targets | Passed across tracked Markdown and new candidate documents |
| Baseline comparison of three demo module sources | Unchanged |
| Development CLI deterministic baseline | Two fresh complete JSON outputs identical; approved fields and numerical versions matched |
| Clean runtime-only wheel install and `pip check` | Passed |
| Installed-wheel demo and module boundary | Passed outside repository root; deterministic baseline matched; only three Python modules installed |
| Linux synthetic container build and execution | Passed; two complete outputs identical and approved non-environment fields matched |
| Git whitespace/diff check | Passed |

The 65 Windows skips comprise 62 Linux filesystem cases, two extended-precision
cases and one unavailable directory-symlink privilege. Two warnings concern the
unchanged historical explicit logistic penalty; no warning-driven scientific
changes were made. Linux full-suite evidence is the successful baseline CI job,
not a newly executed local full Linux suite. Only documentation changed afterward.
Validated profiles remain glibc Linux x86_64 / Python 3.11.17 and Windows AMD64 /
Python 3.11.9. macOS, ARM and musl locks are not supported.

The separate clean Windows runtime installed 12 hash-locked dependencies and the
built wheel with `--no-deps`. Checks ran outside the repository root and verified
that imports resolve inside that environment, that its only package Python files
are the three allowlisted modules, and that research evaluation/modeling modules
are absent. The default demo retained seed 20261026, 100/40 synthetic samples,
12 invented genes and five classes. Its synthetic accuracy/macro F1 of 1.0 are
teaching-fixture outputs, not the held-out scientific figures.

The Linux image was built locally from the unchanged digest-pinned Dockerfile.
It runs as UID/GID 65532:65532. Two runs used `--network none --read-only
--cap-drop ALL --security-opt no-new-privileges --tmpfs /tmp:rw,nosuid,size=64m`,
with no host volumes. Their complete output matched and approved baseline fields
matched; runtime Python was 3.11.17. No image was uploaded.

Initial checks in the pre-existing environment found the missing locked `build`
dependency (334 passed, one dependency-closure failure); that environment was
replaced for validation, without changing locks. Windows' managed default temporary
directory also denied metadata-file renames; workspace-local TEMP/TMP resolved
the editable install and build, as already documented in Phase 4D. These initial
environment failures are not counted as successful candidate checks. Hash-locked
downloads used system certificate trust with TLS verification retained.

Final candidate artifacts (local `dist/`, intentionally ignored by Git):

| Archive | SHA-256 |
|---|---|
| `oncorna_ml-0.1.0-py3-none-any.whl` | `2a5b6cab6d5662809e1b53dedf5779e9e2b9569623124223c289c0d0d982418a` |
| `oncorna_ml-0.1.0.tar.gz` | `4e4762f50ade7bbb1349c9c040bf94f3778cf477012e304e3f85cb390ca23475` |

Metadata inspection confirmed a pure-Python `py3-none-any` wheel, setuptools
84.0.0, the three documented runtime requirements and license file. The archive
tag does not extend support beyond the two locked native dependency profiles.
The source archive intentionally omits full repository documentation; README links
to historical reports and candidate notes are resolved from the repository.
Audits reject unexpected members, unsafe paths, links, duplicate entries, binary
or oversized payloads, recognizable credentials and private paths. Extraction
copied only already-audited regular files into checked destinations.

Recheck from the locked development environment with `python -m pip check`,
`python -m ruff check .`, `python -m ruff format --check src scripts tests tools`,
`python -B -m pytest -q -p no:cacheprovider --basetemp=.pytest_cache/phase5-tests`,
`python -B tools/check_demo.py`, and the build/audit commands in
[packaging](../docs/packaging.md). Clean-wheel installation and restricted-container
commands are in [candidate notes](../docs/release_candidate.md) and that guide.
Do not rebuild these recorded artifacts casually: archive bytes may differ across
builds even when source and scientific behavior are unchanged.

## Review decision and subsequent release

GO for independent final release-candidate review after the checks recorded above.
NO-GO for publication until independent review, exact-final-commit Linux/Windows CI,
third-party preview redistribution review and separate release authorization.
The local commit SHA is supplied in the handoff (a commit cannot contain its own SHA).
The working tree must be clean after the local commit. No push, tag, GitHub release,
package upload or container-image publication is performed.

The exact subsequently authorized release sequence is in
[candidate notes](../docs/release_candidate.md#steps-after-separate-release-authorization).
Adding the actual citation date will require rebuilding and recording replacement
hashes; these candidate hashes must not be presented as checksums of changed archives.
