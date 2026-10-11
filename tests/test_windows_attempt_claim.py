"""Synthetic-only checks of the Windows local claim primitive."""

from __future__ import annotations

import os
from pathlib import Path

import pytest

import oncorna.windows_attempt_claim as claim_module
from oncorna.windows_attempt_claim import create_and_verify_claim


@pytest.mark.skipif(os.name != "nt", reason="Windows claim primitive")
def test_windows_claim_exclusive_content_and_acl(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    path = tmp_path / "synthetic_attempt.json"
    payload = b'{"synthetic":true}\n'
    original = claim_module._file_identity
    observed: list[tuple[int, bytes]] = []

    def record_identity(kernel: object, handle: int) -> tuple[int, bytes]:
        identity = original(kernel, handle)
        observed.append(identity)
        return identity

    monkeypatch.setattr(claim_module, "_file_identity", record_identity)
    create_and_verify_claim(path, payload)
    assert len(observed) == 2  # Creating handle, then independently opened path.
    assert observed[0] == observed[1]
    assert path.read_bytes() == payload
    with pytest.raises(FileExistsError):
        create_and_verify_claim(path, payload)
    assert path.read_bytes() == payload


@pytest.mark.skipif(os.name != "nt", reason="Windows claim primitive")
def test_windows_claim_rejects_substituted_path(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    path = tmp_path / "synthetic_replaced_attempt.json"
    monkeypatch.setattr(claim_module, "_path_identity", lambda *_: (0, b"changed"))
    with pytest.raises(ValueError, match="path identity changed"):
        create_and_verify_claim(path, b"synthetic\n")
    assert path.exists()
    with pytest.raises(FileExistsError):
        create_and_verify_claim(path, b"synthetic\n")


@pytest.mark.skipif(os.name != "nt", reason="Windows claim primitive")
def test_windows_identity_query_failure_closes_handle_and_frees_descriptors(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    path = tmp_path / "synthetic_identity_failure.json"
    original_close = claim_module._close_handle
    original_free = claim_module._free_descriptor
    closed: list[int] = []
    freed: list[object] = []

    def close(kernel: object, handle: int) -> None:
        closed.append(handle)
        original_close(kernel, handle)

    def free(kernel: object, descriptor: object) -> None:
        freed.append(descriptor)
        original_free(kernel, descriptor)

    monkeypatch.setattr(claim_module, "_close_handle", close)
    monkeypatch.setattr(claim_module, "_free_descriptor", free)
    monkeypatch.setattr(
        claim_module,
        "_file_identity",
        lambda *_: (_ for _ in ()).throw(OSError("synthetic identity query failure")),
    )
    with pytest.raises(OSError, match="identity query failure"):
        create_and_verify_claim(path, b"synthetic\n")
    assert len(closed) == 1
    assert len(freed) == 4
    assert path.read_bytes() == b"synthetic\n"
    with pytest.raises(FileExistsError):
        create_and_verify_claim(path, b"synthetic\n")


@pytest.mark.skipif(os.name != "nt", reason="Windows claim primitive")
def test_windows_path_identity_query_failure_retains_claim(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    path = tmp_path / "synthetic_path_query_failure.json"
    monkeypatch.setattr(
        claim_module,
        "_path_identity",
        lambda *_: (_ for _ in ()).throw(OSError("synthetic path identity failure")),
    )
    with pytest.raises(OSError, match="path identity failure"):
        create_and_verify_claim(path, b"synthetic\n")
    assert path.read_bytes() == b"synthetic\n"
    with pytest.raises(FileExistsError):
        create_and_verify_claim(path, b"synthetic\n")


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
