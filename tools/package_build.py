"""Setuptools command limiting installed code to the reviewed synthetic interface."""

from __future__ import annotations

import json
from pathlib import Path

from setuptools.command.build_py import build_py


class SafeBuildPy(build_py):
    def _modules(self) -> set[str]:
        root = Path(__file__).resolve().parents[1]
        return set(json.loads((root / "packaging/distribution-policy.json").read_text())["modules"])

    def find_package_modules(self, package: str, package_dir: str) -> list:
        if package != "oncorna":
            raise ValueError("Only the reviewed oncorna demo package may be built")
        allowed = self._modules()
        return [
            entry
            for entry in super().find_package_modules(package, package_dir)
            if Path(entry[2]).name in allowed
        ]

    def run(self) -> None:
        allowed = {f"oncorna/{name}" for name in self._modules()}
        build_root = Path(self.build_lib)
        if build_root.exists():
            for path in build_root.rglob("*"):
                if path.is_file() and (
                    path.is_symlink() or path.relative_to(build_root).as_posix() not in allowed
                ):
                    raise ValueError(
                        "Unapproved stale build output; use a fresh code-only build tree"
                    )
        super().run()
