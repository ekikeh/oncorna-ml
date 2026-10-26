# Phase 5C — final provenance and v0.1.0 release controls

Preparation date: 2026-10-27. Proposed tag: `v0.1.0`; package/citation version:
`0.1.0`. This checkpoint is local preparation, not a release or publication.
The final local commit SHA is supplied in the handoff; no commit can contain its own SHA.

## Starting state and scope

Local HEAD and live GitHub main were both
`6ccfc76156933c7b5564547c762f25994e79df16`, with a clean working tree.
Read-only GitHub API checks confirmed no `v0.1.0` tag or release (both endpoints
returned 404). [Hosted run 40390676509](https://github.com/ekikeh/oncorna-ml/actions/runs/40390676509)
completed successfully on that baseline for Ubuntu 24.04 / Python 3.11.17 and
Windows 2025 / Python 3.11.9. Its previously inspected logs report Linux 387 passed,
13 skipped, 27 warnings and Windows 336 passed, 64 skipped, 2 warnings. That run
does not cover this later preparation commit.

Exactly five documentation files change:

1. `docs/data_provenance.md`: sourced redistribution evidence, qualifications,
   attribution/privacy conditions and preservation of the historical retrieval date.
2. `docs/demo_dataset.md`: replace the unresolved-license caveat with a link to
   the evidence and its boundaries.
3. `data/demo/README.md`: align preview attribution, privacy and distribution scope.
4. `docs/release_candidate.md`: record baseline CI, provenance resolution, actual-date
   policy and exact-commit/artifact/snapshot authorization controls.
5. `reports/phase_5c_release_controls.md`: this evidence and final gate report.

No scientific source, frozen result, preview values/manifest, synthetic generator,
prediction contract, dependency lock, package metadata, CI, Docker or distribution
policy changes. CITATION.cff remains unchanged, with author handle `ekikeh`, version
0.1.0 and no `date-released`. No new research computation or protected-data access
occurred. The held-out evaluation remains closed.

## Redistribution and privacy decision

All four requested primary sources were reachable and independently checked:

- [Expression metadata](https://tcga.xenahubs.net/download/TCGA.BRCA.sampleMap/HiSeqV2.json):
  resource `TCGA.BRCA.sampleMap/HiSeqV2`, version 2017-10-13, TCGA owner,
  public security and explicit `redistribution: true`; UNC producer and UCSC
  processing are identified. This field is specific to the expression resource.
- [NCI correspondence](https://github.com/cognoma/cancer-data/files/370246/TCGA-public-domain-email.pdf):
  the 2016-07-18 statement by Amy Blum permits redistribution of public-domain,
  open-access TCGA data and requests TCGA attribution. The PDF is a recipient-hosted
  historical correspondence record, not a new formal dataset license; controlled
  data and privacy obligations are outside its permission claim.
- [Xena staff clarification](https://github.com/cognoma/cancer-data/issues/3#issuecomment-257434747):
  the 2016-10-31 statement by Jing Zhu imposes no added restrictions on processed
  TCGA data and treats it as public domain; it explicitly distinguishes other providers.
- [GDC access policy](https://gdc.cancer.gov/access-data/data-access-policies):
  NIH sharing-policy compliance includes no participant identification and attribution
  of datasets/accessions and source repositories. Controlled data have separate rules.

Project interpretation: these statements support redistribution of this documented
open-TCGA training preview, including source-derived labels. The clinical metadata
does not itself expose `redistribution: true`; the broader TCGA/Xena statements and
documented source lineage support that portion of the assessment. Public availability
or attribution alone is not treated as permission. The evidence does not grant rights
over unrelated Xena resources, relicense data under MIT/CC0, certify anonymity or
authorize clinical use. No contradictory statement was found in the reviewed sources.

Retain TCGA, UCSC Xena and UNC attribution, dataset names/versions/URLs, original
checksum references and deterministic preview-selection/value semantics. Identifiers
are linkable, not anonymous; do not attempt reidentification. No participant-level
values are reproduced here. The [provenance document](../docs/data_provenance.md)
distinguishes the explicit statements from the interpretation and links the evidence.

## Complete repository snapshot review

The full baseline tree contains 141 regular files. This preparation adds only this
report, giving 142 intended files. The complete tree and actual in-memory `git archive`
inventory were inspected; there were no links or special files. After the local
commit, the exact inventory is reproducible with `git ls-tree -r --name-only HEAD`.

| Area | Final files | Reviewed contents |
|---|---:|---|
| Repository root | 15 | README, citation/license, packaging, build and editor controls |
| `.github/` | 3 | Issue templates and synthetic CI workflow |
| `configs/` | 6 | Historical scientific/demo settings, not execution authorizations |
| `data/` | 7 | Source checksum reference, four preview files, two empty placeholders |
| `docs/` | 16 | Methodology, scope, provenance, model card and demo/release guidance |
| `packaging/` | 2 | Exact distribution policy and approved synthetic baseline |
| `reports/` | 27 | Historical/aggregate reports, this report and four existing QC figures |
| `requirements/` | 7 | Reviewed input constraints and native platform hash locks |
| `scripts/` | 12 | Historical research/maintenance source |
| `src/oncorna/` | 19 | Historical research source and the three installable demo modules |
| `tests/` | 22 | Synthetic/mocked and platform test infrastructure |
| `tools/` | 5 | Package build, audits, lock and demo checks |
| `workflow/` | 1 | Future workflow plan |

The exact `data/` inventory is `MANIFEST.sha256`, `demo/README.md`,
`demo/demo_expression.tsv`, `demo/demo_metadata.tsv`, `demo/demo_manifest.json`,
`raw/.gitkeep` and `processed/.gitkeep`. Only the documented training preview is
included; raw/processed placeholders are empty. Preview payloads and existing PNG
figures were not parsed or re-derived; baseline blob identity was preserved. The
full tracked text outside preview payloads was scanned for recognizable credentials
and private workspace paths without findings. Historical source/report references
to controls are distinct from real authorization, claim or sealed artifact files.
No fitted binaries, full patient arrays, held-out prediction files, credentials,
private authorization records, claims, seals, environments or caches occur in the
reviewed snapshot inventory. Pattern scans do not prove absence of every possible secret.

Verdict: GO for the reviewed repository snapshot, with the preview attribution and
non-reidentification conditions. GitHub-generated snapshots include the preview and
historical source/tests. They are not the 31-file installable source archive.
The eventual exact release snapshot must be checked before tagging and again after
GitHub generates it; no version tag or snapshot publication was created here.

## Candidate build, installation and checks

All 28 explicitly allowlisted source inputs were copied from immutable baseline Git
blobs into a new clean build tree and compared with current candidate inputs.
They are unchanged by the five documentation edits, which are outside the package
allowlist. This avoids stale build output and excludes all real data from the build.
Native CPython 3.11.9 / Windows AMD64 used the existing verified development lock;
`pip check` passed. `python -m build --no-isolation` built the source archive and
then the wheel from that archive. Workspace-local TEMP/TMP were used as previously
documented. No dependency or backend version changed.

The exact archive inventory/content auditor passed: eight wheel files and 31 source
files. Version/name, Python constraint, MIT metadata, expected three modules and
absence of research modules were verified. Only audited regular payloads were copied
to dedicated resolved extraction destinations beneath their extraction roots.
Traversal, duplicate, link, oversized/binary and recognizable sensitive-content
checks passed. The extracted source archive's 30 focused tests passed.

Newly built artifacts are preserved locally in `dist/phase5c/` with `SHA256SUMS.txt`:

| Archive | SHA-256 |
|---|---|
| `oncorna_ml-0.1.0-py3-none-any.whl` | `0087fdebca1a07ea893e14598b1d7dc858363bef9a037fb00fd02859eeb0465c` |
| `oncorna_ml-0.1.0.tar.gz` | `c75b4816401ef17ec4afcbe3c3b1331707c9285b16ae32d5709fcd0887b9f29c` |

These differ from Phase 5A's wheel `2a5b6cab6d5662809e1b53dedf5779e9e2b9569623124223c289c0d0d982418a`
and source archive `4e4762f50ade7bbb1349c9c040bf94f3778cf477012e304e3f85cb390ca23475`.
The clean source used Git's canonical LF rather than working-copy CRLF in some
files; archive timestamps/compression metadata also differ. Wheel payload differences
are `__init__.py`, LICENSE and their RECORD hashes. Source payload differences are
CITATION.cff, LICENSE, MANIFEST.in, README.md, synthetic-baseline.json, pyproject.toml
and `__init__.py`. Every non-RECORD payload matches Phase 5A after CRLF/LF
normalization; no scientific, executable-logic or version difference was found.
Byte-identical rebuild reproducibility is not claimed. Only the actual attached
artifact bytes and their final hashes can identify an eventual published distribution.

A new dedicated runtime environment installed all 12 runtime-lock packages with
`--require-hashes --no-deps`, then this exact wheel with `--no-deps`.
`pip check` passed. Outside the repository root, two fresh complete demo outputs
matched each other and the approved non-environment baseline. Numerical versions
matched the baseline. Scope remains `synthetic_demo_only`, seed 20261026,
100/40 synthetic samples, 12 invented genes and five classes. Installed Python
files are only `__init__.py`, `prediction.py` and `synthetic_demo.py`; real evaluation
modules are absent. No dataset or fitted research model is needed or loaded.

| Check | Result |
|---|---|
| Focused distribution, lock, prediction and synthetic tests | 137 passed, 2 platform skips, no warnings |
| Extracted source focused tests | 30 passed |
| Ruff and formatting | Passed; all 58 Python files formatted |
| Local documentation targets and new section anchors | Passed |
| Four requested upstream provenance references | Retrieved and content checked |
| Citation, TOML and version consistency | Passed; 0.1.0 and no release date |
| Archive inventories, contents and safe extraction | Passed: 8 / 31 files |
| Fresh runtime wheel smoke/determinism and dependency check | Passed |
| Git diff, scope and whitespace checks | Passed; documentation only |

The two focused skips are Windows `longdouble` exponent-range limitations. No
unchanged full scientific suite or Linux container build was repeated in Phase 5C;
baseline hosted CI provides those results, not CI on this local preparation commit.
Supported locks remain glibc Linux x86_64 / Python 3.11.17 and Windows AMD64 /
Python 3.11.9. macOS, ARM and musl remain unsupported. Historical explicit-logistic-
penalty deprecation warnings remain documented; the new locks do not reproduce
the historical scientific runtime bit for bit.

## Scientific and release gates

Scientific source and approved aggregate report/model card are unchanged. The
task remains existing expression-derived TCGA-BRCA PAM50 reconstruction, not
biomarker discovery or validated clinical prediction. Four-fold development CV is
not nested CV. Held-out macro F1 0.7934, balanced accuracy 0.7556, accuracy 0.8521
(144/169) and Normal-like recall 2/5 remain the approved results. Synthetic metrics
are teaching outputs. External validation, calibration and clinical validity remain
unestablished; historical preprocessing isolation limitations remain disclosed.

GO for independent review and, after separate authorization, publication of this
preparation commit. No unresolved source-access or redistribution contradiction
was found. NO-GO for creating `v0.1.0` now: final hosted CI, date handling, exact
release artifacts and explicit release authorization remain outstanding.

Keep `date-released` omitted until an actual publication date is established and
editing is authorized. Set that real date, review/commit it, then obtain successful
hosted Linux and Windows CI on the exact final publication SHA. If the date changes,
correct the citation and repeat affected checks. Build and audit final archives
from that clean commit, run the clean-wheel/determinism checks, record hashes for
the exact bytes, and inspect the full repository snapshot before tagging.
Obtain explicit approval naming the final SHA, `v0.1.0` and attachments. Only after
that authorization may the tag be created/pushed at that SHA and the GitHub release
created with the audited wheel, package source archive and matching checksum manifest.
Verify uploaded hashes and generated repository snapshots. Package-index or registry
publication requires separate authorization. See the precise
[release procedure](../docs/release_candidate.md#steps-after-separate-release-authorization).

No push, tag, GitHub release, package upload or container publication occurred here.
