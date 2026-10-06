"""Download public TCGA-BRCA expression and clinical data from UCSC Xena.

Run from any working directory with:
    python scripts/download_data.py

Use --force to replace previously downloaded files.
"""

from __future__ import annotations

import argparse
import hashlib
import sys
from dataclasses import dataclass
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

PROJECT_ROOT = Path(__file__).resolve().parents[1]
RAW_DIR = PROJECT_ROOT / "data" / "raw"
MANIFEST_PATH = PROJECT_ROOT / "data" / "MANIFEST.sha256"
CHUNK_SIZE = 1024 * 1024  # 1 MiB; downloads are streamed in chunks.
TIMEOUT_SECONDS = 60
USER_AGENT = "OncoRNA-ML/0.1 (public research data download; https://github.com/ekikeh/oncorna-ml)"


@dataclass(frozen=True)
class SourceFile:
    """One public file to retrieve from UCSC Xena."""

    filename: str
    url: str
    description: str


@dataclass(frozen=True)
class DownloadResult:
    """Basic metadata for a downloaded or already-present file."""

    path: Path
    size_bytes: int
    sha256: str
    downloaded: bool


class DownloadError(RuntimeError):
    """Raised when a source file cannot be downloaded safely."""


# These paired files come from the TCGA Breast Cancer cohort on the UCSC Xena
# TCGA hub. The expression file is the legacy IlluminaHiSeq RNASeqV2 matrix;
# it contains gene-level log2(normalized count + 1) values, not raw integer counts.
SOURCE_FILES = (
    SourceFile(
        filename="TCGA-BRCA_HiSeqV2.tsv.gz",
        url="https://tcga.xenahubs.net/download/TCGA.BRCA.sampleMap/HiSeqV2.gz",
        description="Gene-level RNA-seq expression (Xena HiSeqV2; log2 normalized counts).",
    ),
    SourceFile(
        filename="TCGA-BRCA_BRCA_clinicalMatrix.tsv",
        url="https://tcga.xenahubs.net/download/TCGA.BRCA.sampleMap/BRCA_clinicalMatrix",
        description="Clinical/sample annotations, including PAM50Call_RNAseq where available.",
    ),
)

MANIFEST_HEADER = (
    "# SHA-256 checksums for OncoRNA-ML downloaded inputs.\n"
    "# Paths are relative to the repo root; source URLs are in docs/data_provenance.md.\n"
)


def sha256_file(path: Path, chunk_size: int = CHUNK_SIZE) -> str:
    """Return the SHA-256 checksum of a file without reading it all into memory."""
    digest = hashlib.sha256()
    with Path(path).open("rb") as file_handle:
        while True:
            chunk = file_handle.read(chunk_size)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def _format_size(size_bytes: int) -> str:
    """Format a byte count in both exact bytes and MiB."""
    mebibytes = size_bytes / (1024 * 1024)
    return f"{size_bytes:,} bytes ({mebibytes:.2f} MiB)"


def download_file(
    source: SourceFile,
    destination: Path,
    *,
    force: bool = False,
    chunk_size: int = CHUNK_SIZE,
    timeout: int = TIMEOUT_SECONDS,
) -> DownloadResult:
    """Download one source to disk, streaming bytes and replacing atomically.

    Existing non-empty files are skipped unless ``force`` is true. A temporary
    ``.part`` file prevents an interrupted download from looking complete.
    """
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)

    if destination.exists() and not force:
        if not destination.is_file():
            raise DownloadError(f"Expected a file but found a directory: {destination}")
        size_bytes = destination.stat().st_size
        if size_bytes == 0:
            raise DownloadError(
                f"Existing file is empty: {destination}. Re-download it with --force."
            )
        return DownloadResult(
            path=destination,
            size_bytes=size_bytes,
            sha256=sha256_file(destination, chunk_size=chunk_size),
            downloaded=False,
        )

    if destination.exists() and destination.is_dir():
        raise DownloadError(f"Expected a file but found a directory: {destination}")

    partial_path = destination.with_name(destination.name + ".part")
    partial_path.unlink(missing_ok=True)
    request = Request(source.url, headers={"User-Agent": USER_AGENT})
    digest = hashlib.sha256()
    downloaded_bytes = 0

    try:
        with urlopen(request, timeout=timeout) as response:
            status = getattr(response, "status", 200)
            if status < 200 or status >= 300:
                raise DownloadError(
                    f"Unexpected HTTP status {status} while downloading {source.url}"
                )

            content_length = response.headers.get("Content-Length")
            expected_bytes = int(content_length) if content_length else None

            with partial_path.open("wb") as output:
                while True:
                    chunk = response.read(chunk_size)
                    if not chunk:
                        break
                    output.write(chunk)
                    digest.update(chunk)
                    downloaded_bytes += len(chunk)

        if downloaded_bytes == 0:
            raise DownloadError(f"The server returned an empty file for {source.url}")
        if expected_bytes is not None and downloaded_bytes != expected_bytes:
            raise DownloadError(
                f"Incomplete download for {source.filename}: expected "
                f"{expected_bytes:,} bytes, received {downloaded_bytes:,} bytes"
            )

        partial_path.replace(destination)
    except DownloadError:
        partial_path.unlink(missing_ok=True)
        raise
    except (HTTPError, URLError, TimeoutError, OSError, ValueError) as error:
        partial_path.unlink(missing_ok=True)
        raise DownloadError(
            f"Could not download {source.filename} from {source.url}: {error}"
        ) from error

    return DownloadResult(
        path=destination,
        size_bytes=downloaded_bytes,
        sha256=digest.hexdigest(),
        downloaded=True,
    )


def _read_manifest(manifest_path: Path) -> dict[str, str]:
    """Read existing checksum entries so new downloads do not erase other rows."""
    entries: dict[str, str] = {}
    if not manifest_path.exists():
        return entries

    for line_number, line in enumerate(manifest_path.read_text(encoding="utf-8").splitlines(), 1):
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        parts = stripped.split(maxsplit=1)
        if len(parts) != 2:
            raise ValueError(f"Malformed manifest entry on line {line_number}: {line!r}")
        checksum, relative_path = parts
        if len(checksum) != 64 or any(char not in "0123456789abcdefABCDEF" for char in checksum):
            raise ValueError(f"Invalid SHA-256 checksum on manifest line {line_number}")
        entries[relative_path.strip()] = checksum.lower()
    return entries


def update_manifest(manifest_path: Path, new_entries: dict[str, str]) -> None:
    """Merge checksum entries into a standard, sorted SHA-256 manifest."""
    manifest_path = Path(manifest_path)
    entries = _read_manifest(manifest_path)
    entries.update(new_entries)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)

    temporary_path = manifest_path.with_name(manifest_path.name + ".tmp")
    try:
        with temporary_path.open("w", encoding="utf-8", newline="\n") as output:
            output.write(MANIFEST_HEADER)
            for relative_path, checksum in sorted(entries.items()):
                output.write(f"{checksum}  {relative_path}\n")
        temporary_path.replace(manifest_path)
    finally:
        temporary_path.unlink(missing_ok=True)


def download_files(
    *,
    raw_dir: Path = RAW_DIR,
    manifest_path: Path = MANIFEST_PATH,
    project_root: Path = PROJECT_ROOT,
    force: bool = False,
    source_files: tuple[SourceFile, ...] = SOURCE_FILES,
) -> list[DownloadResult]:
    """Download all configured files and update the repository checksum manifest."""
    raw_dir = Path(raw_dir)
    manifest_path = Path(manifest_path)
    project_root = Path(project_root)
    raw_dir.mkdir(parents=True, exist_ok=True)

    if not raw_dir.is_dir():
        raise DownloadError(f"Could not create the raw-data directory: {raw_dir}")

    results: list[DownloadResult] = []
    for source in source_files:
        destination = raw_dir / source.filename
        result = download_file(destination=destination, source=source, force=force)
        try:
            relative_path = result.path.resolve().relative_to(project_root.resolve()).as_posix()
        except ValueError as error:
            raise ValueError(
                f"Downloaded file {result.path} is outside project root {project_root}"
            ) from error

        # Update after each file so a later failure does not lose earlier checksums.
        update_manifest(manifest_path, {relative_path: result.sha256})
        action = "Downloaded" if result.downloaded else "Already exists; skipped download"
        print(f"{action}: {relative_path}")
        print(f"  description: {source.description}")
        print(f"  size: {_format_size(result.size_bytes)}")
        print(f"  SHA-256: {result.sha256}")
        results.append(result)

    try:
        manifest_display = manifest_path.resolve().relative_to(project_root.resolve()).as_posix()
    except ValueError:
        manifest_display = str(manifest_path)
    print(f"Updated checksum manifest: {manifest_display}")
    return results


def build_parser() -> argparse.ArgumentParser:
    """Create the command-line argument parser."""
    parser = argparse.ArgumentParser(
        description="Download public TCGA-BRCA expression and clinical/PAM50 files from UCSC Xena."
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Download again and replace local files, even if they already exist.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    """Run the downloader and return a shell-friendly exit code."""
    args = build_parser().parse_args(argv)
    print("Source: public UCSC Xena TCGA hub (no credentials or controlled-access files).")

    try:
        download_files(force=args.force)
    except (DownloadError, OSError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
