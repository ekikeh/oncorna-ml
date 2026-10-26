"""Native platform locks must pin and hash the complete installed dependency closure."""

from __future__ import annotations

import importlib.metadata
import re
import sys
from pathlib import Path

import pytest
from packaging.requirements import Requirement
from packaging.utils import canonicalize_name

ROOT = Path(__file__).resolve().parents[1]


def records(path: Path) -> dict[str, str]:
    result = {}
    for line in path.read_text().splitlines():
        if line.startswith(("#", "--")) or not line:
            continue
        match = re.fullmatch(r"([a-z0-9-]+)==([^ ]+) --hash=sha256:([0-9a-f]{64})", line)
        assert match, line
        assert match[1] not in result
        result[match[1]] = match[2]
    return result


def test_platform_locks_are_exact_and_runtime_is_consistent() -> None:
    for target in ("linux", "windows"):
        runtime = records(ROOT / f"requirements/runtime-{target}-py311.lock")
        dev = records(ROOT / f"requirements/dev-{target}-py311.lock")
        assert len(runtime) == 12 and dev
        assert all(dev[name] == version for name, version in runtime.items())
        assert runtime["numpy"] == "2.4.6"
        assert runtime["pandas"] == "2.3.3"
        assert runtime["scikit-learn"] == "1.9.1"
        assert dev["setuptools"] == "84.0.0"


def test_native_lock_matches_installed_versions_and_transitive_requirements() -> None:
    target = "windows" if sys.platform == "win32" else "linux"
    for group in ("runtime", "dev"):
        locked = records(ROOT / f"requirements/{group}-{target}-py311.lock")
        for name, version in locked.items():
            assert importlib.metadata.version(name) == version, name
            for requirement in importlib.metadata.requires(name) or []:
                item = Requirement(requirement)
                if item.marker and not item.marker.evaluate({"extra": ""}):
                    continue
                key = canonicalize_name(item.name)
                assert key in locked, (group, name, requirement)
                assert locked[key] in item.specifier, (name, requirement)


@pytest.mark.parametrize("case", ["platform", "patch", "architecture", "word-size"])
def test_lock_generation_refuses_wrong_profile_before_network(monkeypatch, case: str) -> None:
    from tools import lock_dependencies as locks

    def forbidden_resolution(*_args, **_kwargs):
        pytest.fail("Unsupported platform reached dependency resolution")

    with monkeypatch.context() as guard:
        guard.setattr(locks.sys, "platform", "linux")
        guard.setattr(locks.sys, "version_info", (3, 11, 17))
        guard.setattr(locks.sys, "maxsize", 2**63 - 1)
        guard.setattr(locks.platform, "machine", lambda: "x86_64")
        guard.setattr(locks.subprocess, "run", forbidden_resolution)
        target = "linux"
        if case == "platform":
            target = "windows"
        elif case == "patch":
            guard.setattr(locks.sys, "version_info", (3, 11, 9))
        elif case == "architecture":
            guard.setattr(locks.platform, "machine", lambda: "aarch64")
        else:
            guard.setattr(locks.sys, "maxsize", 2**31 - 1)
        with pytest.raises(ValueError):
            locks.generate(target)
