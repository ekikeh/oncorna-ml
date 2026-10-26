"""Fail-closed wheel/sdist inventory and content inspection without extraction.

Only reviewed text paths are allowed. Inspection never executes package contents.
"""

from __future__ import annotations

import argparse
import json
import re
import stat
import tarfile
import zipfile
from pathlib import Path, PurePosixPath

ROOT = Path(__file__).resolve().parents[1]
POLICY = ROOT / "packaging/distribution-policy.json"
VERSION = "0.1.0"
DIST_INFO = f"oncorna_ml-{VERSION}.dist-info"
EGG_FILES = {"SOURCES.txt"}
WHEEL_FILES = {"METADATA", "WHEEL", "RECORD", "top_level.txt", "licenses/LICENSE"}
FORBIDDEN_PATH = re.compile(
    r"(?:^|/)(?:data|\.git|\.env|\.aws|\.ssh|__pycache__|\.venv|sealed|claims)(?:/|$)"
    r"|\.(?:joblib|pkl|pickle|npy|npz|tsv|csv|pem|key|pyc)$",
    re.IGNORECASE,
)
FORBIDDEN_CONTENT = re.compile(
    rb"TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}"
    rb"|-----BEGIN (?:RSA |OPENSSH |EC )?PRIVATE KEY"
    rb"|gh[pousr]_[A-Za-z0-9]{20,}|AKIA[A-Z0-9]{16}"
    rb"|[A-Za-z]:[\\/]Users[\\/][^\s]+|/(?:home|Users)/[^\s/]+/"
    rb'|"attempt_id"\s*:\s*"[0-9a-f]{8}-[0-9a-f-]{27,}"'
    rb"""|(?im:(?:^\s*|[,{]\s*)["']?(?:api_key|password|access_token|client_secret)["']?\s*[:=]\s*["'][^"']+["'])"""
)


class DistributionError(ValueError):
    """An archive violates the reviewed publication boundary."""


def allowed_files(kind: str, policy: dict | None = None) -> set[str]:
    policy = policy if policy is not None else json.loads(POLICY.read_text())
    if kind == "wheel":
        return {f"oncorna/{name}" for name in policy["modules"]} | {
            f"{DIST_INFO}/{name}" for name in WHEEL_FILES
        }
    if kind != "sdist":
        raise DistributionError("Unknown archive kind")
    return (
        set(policy["sdist_files"])
        | {f"src/oncorna/{name}" for name in policy["modules"]}
        | {f"src/oncorna_ml.egg-info/{name}" for name in EGG_FILES}
        | {"PKG-INFO", "setup.cfg"}
    )


def check_member(name: str, content: bytes, allowed: set[str]) -> None:
    path = PurePosixPath(name)
    if path.is_absolute() or ".." in path.parts or "\\" in name or ":" in name:
        raise DistributionError("Unsafe archive path")
    if str(path) != name or FORBIDDEN_PATH.search(name) or name not in allowed:
        raise DistributionError(f"Unapproved archive path: {name}")
    if len(content) > 2_000_000 or b"\x00" in content:
        raise DistributionError(f"Binary or oversized content: {name}")
    try:
        content.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise DistributionError(f"Non-text content: {name}") from exc
    if FORBIDDEN_CONTENT.search(content):
        raise DistributionError(f"Prohibited sensitive content: {name}")


def audit(path: Path) -> list[str]:
    """Require exactly the approved inventory and reject links/duplicates/content leaks."""
    kind = "wheel" if path.suffix == ".whl" else "sdist"
    allowed = allowed_files(kind)
    contents = {}

    def add(name: str, payload: bytes) -> None:
        if name in contents:
            raise DistributionError("Duplicate archive entry")
        check_member(name, payload, allowed)
        contents[name] = payload

    if kind == "wheel":
        with zipfile.ZipFile(path) as archive:
            for member in archive.infolist():
                if member.is_dir():
                    raise DistributionError("Wheel contains unexpected directory entry")
                if member.file_size > 2_000_000 or member.flag_bits & 1:
                    raise DistributionError("Oversized or encrypted archive entry")
                if stat.S_ISLNK(member.external_attr >> 16):
                    raise DistributionError("Archive link is forbidden")
                add(member.filename, archive.read(member))
    else:
        prefix = f"oncorna_ml-{VERSION}/"
        allowed_directories = {
            parent.as_posix()
            for name in allowed
            for parent in PurePosixPath(name).parents
            if parent.as_posix() != "."
        }
        with tarfile.open(path, "r:gz") as archive:
            for member in archive:
                if not member.name.startswith(prefix) and member.name != prefix.rstrip("/"):
                    raise DistributionError("Unexpected source archive root")
                relative = member.name[len(prefix) :]
                if member.isdir():
                    directory = relative.rstrip("/")
                    if directory and (
                        directory not in allowed_directories or FORBIDDEN_PATH.search(directory)
                    ):
                        raise DistributionError("Unapproved archive directory")
                    continue
                if not member.isfile() or member.size > 2_000_000:
                    raise DistributionError("Archive link, special file or oversized entry")
                handle = archive.extractfile(member)
                if handle is None:
                    raise DistributionError("Unreadable source archive member")
                add(relative, handle.read())
    if set(contents) != allowed:
        raise DistributionError(
            f"Incomplete archive inventory: {sorted(allowed - contents.keys())}"
        )
    return sorted(contents)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archives", nargs="+", type=Path)
    args = parser.parse_args()
    for archive in args.archives:
        names = audit(archive)
        print(json.dumps({"archive": archive.name, "file_count": len(names), "files": names}))
