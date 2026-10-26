# Phase 4D — packaging and synthetic-environment validation

Reviewed starting revision: `eeeeae8e268059274dbce4423489c38fc43c9c11`.
This is an engineering checkpoint for the public synthetic demo, not a new scientific
analysis. No protected research input, fitted artifact or final-test result was opened.

## Baseline audit and resulting scope

Before configuration changes, a sanitized code-only snapshot of the starting revision
produced a 24-file wheel (19 modules plus metadata/license files) and a 49-file source
archive. The source archive included tests without their full repository configuration.
The snapshot deliberately omitted all data and the real-training preview; its inventory
does not imply that those assets were packaged or inspected.

Content inspection of the broad library wheel detected patient-like identifiers in
historical source. The public wheel is therefore deliberately limited to the unchanged
`__init__.py`, `prediction.py` and `synthetic_demo.py` modules. Research source, split and
evaluation safeguards remain unchanged in the repository. The resulting wheel has
**8 approved text files** and the source archive **31 approved text files**. Both archives
pass the exact allowlist and forbidden-content audit. Linux archives were safely extracted
after passing audit; the source archive's 30 focused packaging tests also passed.

## Dependency and installation evidence

- Runtime: 12 exact, SHA-256-locked packages on each platform.
- Development: 32 packages on Linux CPython 3.11.17 x86_64; 33 on Windows CPython 3.11.9
  AMD64, including colorama.
- Native pip resolution used public PyPI and reviewed direct/transitive constraints.
- Clean installations and `pip check` passed. Tests verify both runtime and development
  dependency closure against installed metadata, exact versions and markers.
- A separate runtime-only environment reinstalled the audited wheel on each platform.
  Neither installed package contains the research evaluation modules.
- Two fresh complete JSON outputs matched within each tested environment and matched
  approved Phase 4C non-environment fields. Numerical dependency versions did not change.

The default seed remains 20261026; training/evaluation populations remain 100/40 samples,
with 12 invented genes and five synthetic classes. No generator or numeric-validator
source changed in Phase 4D. Cross-environment bitwise equality is not claimed.

## Test and configuration evidence

| Check | Local result |
|---|---|
| Full isolated Linux Python 3.11 suite | 387 passed, 13 Windows-only skips, 27 historical warnings |
| Full code-only Windows Python 3.11 suite | 335 passed, 65 skips, 2 historical warnings |
| Packaging/lock-focused tests | 30 passed |
| Extracted source archive's focused suite | 30 passed |
| Ruff and repository Python formatting | Passed |
| TOML/YAML, local documentation links, lock-input digests | Passed |
| Wheel/source inventories and content audits | Passed |
| Same-environment synthetic scope/numerical/determinism checks | Passed |

Linux skips are the 13 Windows ACL/guard cases covered on Windows. Windows skips comprise
62 Linux filesystem-integration cases, two extended-precision cases, and one directory
symlink privilege limitation. Linux exercises the corresponding filesystem and extended-
precision behavior. Warnings concern the unchanged historical explicit logistic penalty;
there were no new convergence or dependency-conflict failures.

The new GitHub Actions workflow has **not run on GitHub** during this checkpoint. Local
validation is not reported as a CI pass. Actions use verified immutable revisions, exact
Python patches, read-only permissions, no persisted checkout credential, and no secrets,
research downloads, package publication, tags or releases.

## Container and environmental checks

The runtime image was built from the documented digest-pinned Python base. Its effective
Docker context contained exactly six approved files; unrelated synthetic sentinel files
were excluded by the actual Docker context exporter. The image contains only three demo
modules plus locked public dependencies, not the repository's research modules or data.

Two fresh runtime containers were inspected and executed: UID/GID 65532, no network,
read-only root, dropped capabilities, no new privileges, disposable temporary storage,
and no host/production/Docker-socket mounts. Complete JSON repeated exactly and approved
synthetic values matched.

Resolved local environment issues were reported separately from code tests: Windows pip
needed system certificate trust with TLS verification retained; managed Windows default
temporary-directory renames required workspace-local TEMP/TMP; Linux's default noexec
temporary mount prevented loading NumPy native libraries in a newly created temporary
venv. The wheel smoke check passed with dedicated executable temporary validation storage.
The production-style demo runtime flags were not weakened to accommodate that test.

## Publication and release boundaries

The installable archive is a minimal demo package; the full repository remains the source
for historical research implementation and the complete mocked test suite. The real
training preview, research data, fitted binaries, held-out predictions, authorization
records, claims, seals, credentials, caches and private workspace paths are excluded.
Source-command filtering and archive audits supplement the explicit source manifest;
`.gitignore` alone is not the distribution boundary. The pre-existing download test has
one formatting-only adjustment so CI can check all repository Python formatting.

Only the two documented 64-bit profiles were validated. No macOS, ARM or musl lock is
claimed. Content scanners complement human review; arbitrary secrets cannot be ruled out
by patterns alone. Independent Phase 4D publication review and actual successful GitHub
CI are still required. Phase 5 should choose the release version, update version-coupled
archive paths/policy, rebuild and audit exact artifacts, record hashes and obtain explicit
release approval. No push, package publication, version tag or release occurred here.

## Exact implementation inventory

29 changed/new tracked-intended files:

```text
.dockerignore
.github/workflows/ci.yml
.gitignore
Dockerfile
MANIFEST.in
README.md
docs/packaging.md
docs/release_scope.md
docs/reproducibility.md
docs/synthetic_demo.md
packaging/distribution-policy.json
packaging/synthetic-baseline.json
pyproject.toml
reports/phase_4d_packaging_validation.md
requirements/dev-linux-py311.lock
requirements/dev-windows-py311.lock
requirements/dev.in
requirements/runtime-linux-py311.lock
requirements/runtime-windows-py311.lock
requirements/runtime.in
requirements/versions.txt
tests/test_distribution.py
tests/test_download_data.py
tests/test_lockfiles.py
tools/__init__.py
tools/audit_distribution.py
tools/check_demo.py
tools/lock_dependencies.py
tools/package_build.py
```
