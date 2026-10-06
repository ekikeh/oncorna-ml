from __future__ import annotations

import hashlib
from io import BytesIO
from pathlib import Path
from urllib.error import URLError

import pytest

from scripts import download_data


class FakeResponse(BytesIO):
    """Small in-memory response used to test streaming without network access."""

    def __init__(self, content: bytes) -> None:
        super().__init__(content)
        self.status = 200
        self.headers = {"Content-Length": str(len(content))}


@pytest.fixture
def source_file() -> download_data.SourceFile:
    return download_data.SourceFile(
        filename="tiny.tsv",
        url="https://example.test/tiny.tsv",
        description="test fixture",
    )


def test_sha256_file_matches_expected_digest(tmp_path: Path) -> None:
    file_path = tmp_path / "sample.bin"
    file_path.write_bytes(b"OncoRNA-ML checksum test")

    assert download_data.sha256_file(file_path) == hashlib.sha256(
        b"OncoRNA-ML checksum test"
    ).hexdigest()


def test_download_creates_raw_directory_and_manifest(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    source_file: download_data.SourceFile,
) -> None:
    content = b"sample\tvalue\nSYN-BRCA-0001\t1\n"
    monkeypatch.setattr(
        download_data,
        "urlopen",
        lambda request, timeout: FakeResponse(content),
    )
    raw_dir = tmp_path / "data" / "raw"
    manifest_path = tmp_path / "data" / "MANIFEST.sha256"

    results = download_data.download_files(
        raw_dir=raw_dir,
        manifest_path=manifest_path,
        project_root=tmp_path,
        source_files=(source_file,),
    )

    saved_file = raw_dir / source_file.filename
    assert raw_dir.is_dir()
    assert saved_file.read_bytes() == content
    assert results[0].downloaded is True
    assert results[0].sha256 == hashlib.sha256(content).hexdigest()
    assert f"{results[0].sha256}  data/raw/{source_file.filename}" in manifest_path.read_text()
    assert not saved_file.with_name(saved_file.name + ".part").exists()


def test_existing_file_is_skipped_without_network_request(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    source_file: download_data.SourceFile,
) -> None:
    raw_dir = tmp_path / "data" / "raw"
    raw_dir.mkdir(parents=True)
    saved_file = raw_dir / source_file.filename
    saved_file.write_bytes(b"already downloaded")

    def unexpected_request(request, timeout):
        raise AssertionError("urlopen should not be called for an existing file")

    monkeypatch.setattr(download_data, "urlopen", unexpected_request)
    result = download_data.download_file(source_file, saved_file)

    assert result.downloaded is False
    assert result.sha256 == hashlib.sha256(b"already downloaded").hexdigest()


def test_force_replaces_existing_file(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    source_file: download_data.SourceFile,
) -> None:
    raw_dir = tmp_path / "data" / "raw"
    raw_dir.mkdir(parents=True)
    saved_file = raw_dir / source_file.filename
    saved_file.write_bytes(b"old content")
    new_content = b"new content"
    monkeypatch.setattr(
        download_data,
        "urlopen",
        lambda request, timeout: FakeResponse(new_content),
    )

    result = download_data.download_file(source_file, saved_file, force=True)

    assert result.downloaded is True
    assert saved_file.read_bytes() == new_content
    assert result.sha256 == hashlib.sha256(new_content).hexdigest()


def test_failed_download_has_clear_error_and_removes_partial_file(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    source_file: download_data.SourceFile,
) -> None:
    destination = tmp_path / "data" / "raw" / source_file.filename

    def failed_request(request, timeout):
        raise URLError("test network failure")

    monkeypatch.setattr(download_data, "urlopen", failed_request)

    with pytest.raises(download_data.DownloadError, match="Could not download tiny.tsv"):
        download_data.download_file(source_file, destination)

    assert not destination.exists()
    assert not destination.with_name(destination.name + ".part").exists()
