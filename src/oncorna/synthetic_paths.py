"""Fail-closed path boundary for synthetic final-test files."""

from __future__ import annotations

import hashlib
import os
import stat
from pathlib import Path

PROJECT_DATA = Path(__file__).resolve().parents[2] / "data"
FileIdentity = tuple[int, int, int, int, int]


class SyntheticIsolationError(PermissionError):
    """Filesystem-backed synthetic execution requires an isolated Linux container."""


def require_isolated_synthetic_filesystem() -> None:
    """Fail before path inspection on the Windows production host."""
    if os.name == "nt":
        raise SyntheticIsolationError(
            "Synthetic filesystem execution is disabled on Windows; use the code-only container"
        )


def _normal(path: Path) -> str:
    return os.path.normcase(os.path.abspath(path))


def _inside(path: Path, root: Path) -> bool:
    location, boundary = _normal(path), _normal(root)
    try:
        return os.path.commonpath((location, boundary)) == boundary
    except ValueError:  # Different Windows drives.
        return False


def _canonical(path: Path, *, must_exist: bool, protected_data: Path) -> Path:
    require_isolated_synthetic_filesystem()
    try:
        resolved = path.resolve(strict=must_exist)
        boundary = protected_data.resolve(strict=True)
    except (OSError, RuntimeError) as exc:
        raise ValueError("Synthetic path identity cannot be verified") from exc
    if _normal(path) != _normal(resolved):
        raise ValueError("Synthetic path uses a symlink, junction or redirection")
    if _inside(resolved, boundary):
        raise ValueError("Synthetic path enters protected production data")
    return resolved


def synthetic_output_path(path: Path, *, protected_data: Path = PROJECT_DATA) -> Path:
    """Check before mkdir; all production data paths are reserved from synthetic output."""
    resolved = _canonical(path, must_exist=False, protected_data=protected_data)
    if path.exists():
        raise FileExistsError(f"Synthetic output already exists: {path}")
    if not path.parent.is_dir():
        raise ValueError("Synthetic output requires an existing parent")
    return resolved


def _identity(path: Path) -> FileIdentity:
    require_isolated_synthetic_filesystem()
    info = path.stat()
    if not stat.S_ISREG(info.st_mode) or not info.st_ino or info.st_nlink != 1:
        raise ValueError("Synthetic input is not a unique regular file")
    return (info.st_dev, info.st_ino, info.st_size, info.st_mtime_ns, info.st_ctime_ns)


def synthetic_input_path(path: Path, *, protected_data: Path = PROJECT_DATA) -> FileIdentity:
    """Validate location and identity without opening or reading the candidate."""
    resolved = _canonical(path, must_exist=True, protected_data=protected_data)
    return _identity(resolved)


def verify_synthetic_input(
    path: Path, expected: FileIdentity, *, protected_data: Path = PROJECT_DATA
) -> None:
    if synthetic_input_path(path, protected_data=protected_data) != expected:
        raise ValueError("Synthetic input changed after path validation")


def verify_synthetic_handle(handle: object, expected: FileIdentity) -> None:
    info = os.fstat(handle.fileno())  # type: ignore[attr-defined]
    actual = (info.st_dev, info.st_ino, info.st_size, info.st_mtime_ns, info.st_ctime_ns)
    if actual != expected:
        raise ValueError("Synthetic input changed at file open")


def sha256_synthetic_input(
    path: Path, expected: FileIdentity, *, protected_data: Path = PROJECT_DATA
) -> str:
    """Hash only a verified open handle, and reject substitution during the read."""
    verify_synthetic_input(path, expected, protected_data=protected_data)
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        verify_synthetic_handle(handle, expected)
        verify_synthetic_input(path, expected, protected_data=protected_data)
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
        verify_synthetic_handle(handle, expected)
    verify_synthetic_input(path, expected, protected_data=protected_data)
    return digest.hexdigest()
