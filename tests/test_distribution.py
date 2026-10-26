"""Generated archives only: regression checks for the public file/content boundary."""

from __future__ import annotations

import io
import json
import stat
import tarfile
import zipfile
from pathlib import Path

import pytest

from tools.audit_distribution import DistributionError, allowed_files, audit


def wheel(tmp_path: Path, changes: dict[str, bytes] | None = None, missing: bool = False) -> Path:
    files = dict.fromkeys(allowed_files("wheel"), b"safe text\n")
    files.update(changes or {})
    if missing:
        files.pop("oncorna/prediction.py")
    output = tmp_path / "fixture.whl"
    with zipfile.ZipFile(output, "w") as archive:
        for name, data in files.items():
            archive.writestr(name, data)
    return output


def test_exact_safe_wheel_inventory(tmp_path: Path) -> None:
    assert set(audit(wheel(tmp_path))) == allowed_files("wheel")


def test_public_modules_are_only_the_demo_interface() -> None:
    modules = {name for name in allowed_files("wheel") if name.startswith("oncorna/")}
    assert modules == {"oncorna/__init__.py", "oncorna/prediction.py", "oncorna/synthetic_demo.py"}


def test_build_refuses_unapproved_stale_outputs(tmp_path: Path) -> None:
    from setuptools import Distribution

    from tools.package_build import SafeBuildPy

    command = SafeBuildPy(Distribution())
    command.build_lib = str(tmp_path)
    (tmp_path / "unreviewed.py").write_text("synthetic stale output\n")
    with pytest.raises(ValueError, match="stale build output"):
        command.run()


@pytest.mark.parametrize(
    "path",
    [
        "data/demo/demo_expression.tsv",
        "data/processed/model.joblib",
        "/outside.py",
        "../outside.py",
        "oncorna/extra.py",
        "oncorna/.env",
        "oncorna/__pycache__/cached.pyc",
        "C:\\secret.py",
    ],
)
def test_unapproved_paths_fail(tmp_path: Path, path: str) -> None:
    with pytest.raises(DistributionError):
        audit(wheel(tmp_path, {path: b"synthetic test sentinel"}))


@pytest.mark.parametrize(
    "payload",
    [
        b"TC" + b"GA-AA-" + b"1234",
        b"-----BEGIN " + b"PRIVATE KEY-----",
        b"ghp_" + b"A" * 30,
        b"C:" + b"\\Us" + b"ers\\SYN\\private",
        b"api_key = " + b'"synthetic-secret"',
        b"\x00binary",
        b"a" * 2_000_001,
    ],
    ids=["patient", "private-key", "token", "private-path", "credential", "binary", "oversized"],
)
def test_sensitive_binary_and_large_content_fail(tmp_path: Path, payload: bytes) -> None:
    with pytest.raises(DistributionError):
        audit(wheel(tmp_path, {"oncorna/synthetic_demo.py": payload}))


def test_missing_required_file_fails(tmp_path: Path) -> None:
    with pytest.raises(DistributionError, match="Incomplete"):
        audit(wheel(tmp_path, missing=True))


def test_zip_symlink_and_duplicate_fail(tmp_path: Path) -> None:
    path = wheel(tmp_path)
    with zipfile.ZipFile(path, "a") as archive:
        link = zipfile.ZipInfo("oncorna/link.py")
        link.create_system = 3
        link.external_attr = (stat.S_IFLNK | 0o777) << 16
        archive.writestr(link, "outside")
    with pytest.raises(DistributionError, match="link"):
        audit(path)
    path = wheel(tmp_path)
    with zipfile.ZipFile(path, "a") as archive:
        with pytest.warns(UserWarning, match="Duplicate"):
            archive.writestr("oncorna/prediction.py", "duplicate")
    with pytest.raises(DistributionError, match="Duplicate"):
        audit(path)


def test_source_archive_inventory_and_links(tmp_path: Path) -> None:
    path = tmp_path / "fixture.tar.gz"
    prefix = "oncorna_ml-0.1.0/"
    with tarfile.open(path, "w:gz") as archive:
        for name in sorted(allowed_files("sdist")):
            content = b"safe text\n"
            item = tarfile.TarInfo(prefix + name)
            item.size = len(content)
            archive.addfile(item, io.BytesIO(content))
    assert set(audit(path)) == allowed_files("sdist")
    with tarfile.open(path, "w:gz") as archive:
        item = tarfile.TarInfo(prefix + "src/oncorna/prediction.py")
        item.type = tarfile.SYMTYPE
        item.linkname = "/outside"
        archive.addfile(item)
    with pytest.raises(DistributionError, match="link"):
        audit(path)


def test_manifest_matches_reviewed_policy() -> None:
    root = Path(__file__).resolve().parents[1]
    policy = json.loads((root / "packaging/distribution-policy.json").read_text())
    lines = (root / "MANIFEST.in").read_text().splitlines()
    assert lines[0] == "global-exclude *"
    assert set(lines[1:]) == {"include " + name for name in policy["sdist_files"]} | {
        "include src/oncorna/" + name for name in policy["modules"]
    }


def test_source_archive_forbids_even_empty_data_directories(tmp_path: Path) -> None:
    path = tmp_path / "empty-data.tar.gz"
    with tarfile.open(path, "w:gz") as archive:
        for name in sorted(allowed_files("sdist")):
            item = tarfile.TarInfo("oncorna_ml-0.1.0/" + name)
            item.size = 1
            archive.addfile(item, io.BytesIO(b"x"))
        directory = tarfile.TarInfo("oncorna_ml-0.1.0/data")
        directory.type = tarfile.DIRTYPE
        archive.addfile(directory)
    with pytest.raises(DistributionError, match="directory"):
        audit(path)


def test_docker_context_has_closed_directory_exceptions() -> None:
    root = Path(__file__).resolve().parents[1]
    expected = [
        "**",
        "!Dockerfile",
        "!.dockerignore",
        "!requirements/",
        "requirements/*",
        "!requirements/runtime-linux-py311.lock",
        "!src/",
        "src/*",
        "!src/oncorna/",
        "src/oncorna/*",
        "!src/oncorna/__init__.py",
        "!src/oncorna/prediction.py",
        "!src/oncorna/synthetic_demo.py",
    ]
    assert (root / ".dockerignore").read_text().splitlines() == expected
