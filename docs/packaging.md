# Reproducible synthetic environments and distribution

These specifications apply to the public `synthetic_demo_only` software fixture.
They do not reconstruct the historical numerical environment of the fitted TCGA-BRCA
model and do not authorize reopening its closed evaluation. Scientific provenance and
limitations remain in the [reproducibility appendix](reproducibility.md).

## Supported environments and exact locks

| Profile | Python | Platform | Runtime lock | Development/test lock |
|---|---|---|---|---|
| Linux | CPython 3.11.17 | glibc Linux x86_64 | `runtime-linux-py311.lock` | `dev-linux-py311.lock` |
| Windows | CPython 3.11.9 | Windows AMD64 | `runtime-windows-py311.lock` | `dev-windows-py311.lock` |

All lock files are under `requirements/`. Each resolved dependency, including transitive
dependencies, has an exact version and the SHA-256 of the selected native wheel.
The runtime closure contains 12 packages. Development includes runtime plus build,
lint and test tools: 32 packages on Linux and 33 on Windows (including colorama).
These are platform-specific wheel locks, not universal Linux/Windows/ARM specifications.
The Windows patch is the published Windows build available in the official
[setup-python version manifest](https://github.com/actions/python-versions/blob/main/versions-manifest.json);
it is a compatibility-test target, separate from the Linux container runtime.

The validated numerical versions remain NumPy 2.4.6, pandas 2.3.3, scikit-learn 1.9.1,
SciPy 1.17.1, joblib 1.6.0 and threadpoolctl 3.7.0. Exact versions for the entire
closure are in the locks. The broad package metadata describes compatibility;
`environment.yml` remains a legacy convenience specification, not an exact lock.

## Clean installation

Create a new environment with the profile's Python; do not reuse the research environment.
For Linux, from the repository root:

```bash
python3.11 -m venv .venv-demo
.venv-demo/bin/python -m pip --isolated install --require-hashes --no-deps -r requirements/dev-linux-py311.lock
.venv-demo/bin/python -m pip --isolated install --no-deps --no-build-isolation -e .
.venv-demo/bin/python -m pip check
.venv-demo/bin/python -B tools/check_demo.py
.venv-demo/bin/python -B -m pytest -q -p no:cacheprovider --basetemp=.pytest_cache/demo-tests
```

On Windows use `py -3.11 -m venv .venv-demo`, `.venv-demo\Scripts\python.exe`, and
`requirements/dev-windows-py311.lock`. If pip 24 cannot see your organization's trusted
root certificate, add `--use-feature=truststore` before `install`; TLS verification
remains enabled. Never replace verification with trusted-host or insecure download flags.

Use the corresponding runtime lock for a demonstration-only install. Install the built
wheel with `--no-deps` afterward. Hash-locked installation plus `pip check` and the
dependency-closure tests prevent unreviewed transitive resolution.

## Deliberate lock regeneration

`requirements/runtime.in` and `dev.in` specify direct pins; `versions.txt` constrains
every reviewed transitive version. On each exact native profile, run:

```bash
python tools/lock_dependencies.py --platform linux
# On the Windows profile instead:
python tools/lock_dependencies.py --platform windows
```

The generator uses pip's native resolver against public PyPI, refuses dependencies not
listed at their reviewed versions, and records selected wheel hashes and an input digest.
It does not install the project. Generation requires network access; running the demo does
not. Regeneration is a reviewed change: inspect every version/hash difference, reinstall
cleanly, run the suite, and compare synthetic outputs before accepting a new lock.
See [pip's repeatable-install guidance](https://pip.pypa.io/en/stable/topics/repeatable-installs/).

`tools/check_demo.py` compares two fresh complete JSON outputs, then checks all non-environment
fields against `packaging/synthetic-baseline.json`, captured from approved Phase 4C.
It also checks the numerical dependency versions. Python patch metadata intentionally
differs between the two profiles. Do not update the baseline to conceal dependency-related
numerical changes; stop for review. No cross-platform or arbitrary-version bitwise promise
is made by same-environment comparisons.

## Distribution boundary

`packaging/distribution-policy.json` enumerates the three reviewed demo modules and
permitted source-archive text files. `MANIFEST.in` first excludes everything, then
explicitly includes those files. Setuptools package discovery is restricted to `oncorna`,
package data is disabled, and build-backend versions are pinned. A small setuptools
build command filters Python modules to that same list and refuses unapproved stale
build output. Historical scientific modules stay in the repository, outside archives.

The wheel contains only `__init__.py`, `prediction.py`, `synthetic_demo.py`, and standard
metadata/license files. The source archive
also contains selected packaging/demo documentation, requirements, maintenance tools,
and two packaging-focused test files. It is not a complete checkout for the historical
scientific workflow or full test suite; those remain available in the repository.
Full cohort files, the real training preview, model binaries, predictions, authorization
records, claims, seals, environment files, credentials and caches are excluded.

```bash
python -m build --no-isolation
python tools/audit_distribution.py dist/oncorna_ml-0.1.0-py3-none-any.whl dist/oncorna_ml-0.1.0.tar.gz
```

Build only after installing the development lock. The auditor inspects archive contents
without executing or extracting them. It requires the exact permitted inventory and
rejects traversal, links, duplicates, binary/oversized payloads, patient identifiers,
recognizable credentials, protected-record payloads and private workspace paths. Automated
patterns supplement human inspection; they cannot establish that arbitrary source text
contains no secret. New files require an explicit policy review. Do not publish an archive
that fails either inventory or content inspection. See [release scope](release_scope.md).

## Minimal container

The Dockerfile uses the verified Python 3.11.17 Linux base at digest
`sha256:b7a9e5f6a93d6d4a31f129bbd4bdee69efccb6ded0fb230d78aafaf37af98377`.
BuildKit caches only public pip downloads; dependency hashes are checked on installation.
`.dockerignore` allows only the runtime lock, Docker configuration and the three demo
modules. Neither patient data nor the full repository enters the intended build context.

```bash
docker build --platform linux/amd64 -t oncorna-demo .
docker run --rm --network none --read-only --cap-drop ALL \
  --security-opt no-new-privileges --tmpfs /tmp:rw,nosuid,size=64m oncorna-demo
```

The image runs as UID/GID 65532 and needs no volumes, credentials or Docker socket.
The network/read-only/capability restrictions are runtime flags; a Dockerfile alone
cannot impose all of them. Temporary writable storage is disposable. This image runs
only the synthetic CLI, not the research evaluation or a web service.

## CI and release gate

The workflow has Linux and Windows jobs with exact Python patches and hash-locked tools,
read-only repository permissions, immutable action references and disabled checkout
credential persistence. It runs lint, formatting, the full applicable synthetic/mocked
suite, deterministic CLI checks, distribution builds/audits, and a wheel smoke test in
a fresh runtime-only environment. Linux also builds and securely runs the minimal image.
No step downloads TCGA data, loads the real model, publishes packages or creates a release.
Platform skips are printed by pytest. The historical explicit-penalty deprecation warnings
remain historical; dependency locks prevent silent upgrades across that compatibility edge.

An implemented workflow is not evidence of a successful GitHub Actions run. Publication
review and a successful actual CI run on the intended release commit remain Phase 5 gates.
