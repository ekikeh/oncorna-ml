"""Synthetic-only checks of the Windows local claim primitive."""

from __future__ import annotations

import os
from pathlib import Path

import pytest

import oncorna.windows_attempt_claim as claim_module
from oncorna.windows_attempt_claim import create_and_verify_claim


@pytest.mark.skipif(os.name != "nt", reason="Windows claim primitive")
def test_windows_claim_exclusive_content_and_acl(tmp_path: Path) -> None:
    path = tmp_path / "synthetic_attempt.json"
    payload = b'{"synthetic":true}\n'
    create_and_verify_claim(path, payload)
    assert path.read_bytes() == payload
    with pytest.raises(FileExistsError):
        create_and_verify_claim(path, payload)
    assert path.read_bytes() == payload


@pytest.mark.skipif(os.name != "nt", reason="Windows claim primitive")
def test_windows_acl_verification_failure_consumes_visible_claim(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    path = tmp_path / "synthetic_acl_failure.json"
    original = claim_module._dacl_string
    calls = 0

    def changed_dacl(descriptor: object) -> str:
        nonlocal calls
        calls += 1
        return original(descriptor) if calls == 1 else "D:P"

    monkeypatch.setattr(claim_module, "_dacl_string", changed_dacl)
    with pytest.raises(PermissionError, match="ACL differs"):
        create_and_verify_claim(path, b"synthetic\n")
    assert path.exists()
    with pytest.raises(FileExistsError):
        create_and_verify_claim(path, b"synthetic\n")


@pytest.mark.skipif(os.name != "nt", reason="Windows claim primitive")
def test_windows_invalid_security_descriptor_fails_before_creation(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    path = tmp_path / "synthetic_acl_creation_failure.json"
    monkeypatch.setattr(claim_module, "_user_sid", lambda: "invalid-SID")
    with pytest.raises(OSError):
        create_and_verify_claim(path, b"synthetic\n")
    assert not path.exists()
