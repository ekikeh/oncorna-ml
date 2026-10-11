"""Windows-only exclusive local claim with an ACL supplied at file creation."""

from __future__ import annotations

import csv
import ctypes
import io
import os
import re
import subprocess
from ctypes import wintypes
from pathlib import Path


class SecurityAttributes(ctypes.Structure):
    _fields_ = [
        ("nLength", wintypes.DWORD),
        ("lpSecurityDescriptor", ctypes.c_void_p),
        ("bInheritHandle", wintypes.BOOL),
    ]


def _user_sid() -> str:
    result = subprocess.check_output(["whoami", "/user", "/fo", "csv", "/nh"], text=True)
    rows = list(csv.reader(io.StringIO(result)))
    if len(rows) != 1 or len(rows[0]) != 2 or not re.fullmatch(r"S-\d+(?:-\d+)+", rows[0][1]):
        raise PermissionError("Cannot establish the Windows claim owner's SID")
    return rows[0][1]


def _dacl_string(security_descriptor: ctypes.c_void_p) -> str:
    advapi = ctypes.WinDLL("advapi32", use_last_error=True)
    kernel = ctypes.WinDLL("kernel32", use_last_error=True)
    convert = advapi.ConvertSecurityDescriptorToStringSecurityDescriptorW
    convert.argtypes = [
        ctypes.c_void_p,
        wintypes.DWORD,
        wintypes.DWORD,
        ctypes.POINTER(ctypes.c_void_p),
        ctypes.POINTER(wintypes.ULONG),
    ]
    convert.restype = wintypes.BOOL
    output = ctypes.c_void_p()
    if not convert(security_descriptor, 1, 4, ctypes.byref(output), None):
        raise ctypes.WinError(ctypes.get_last_error())
    try:
        return ctypes.wstring_at(output)
    finally:
        kernel.LocalFree.argtypes = [ctypes.c_void_p]
        kernel.LocalFree(output)


def create_and_verify_claim(path: Path, payload: bytes) -> None:
    """Create one claim, flush and verify its bytes and protected DACL on its handle."""
    if os.name != "nt" or not path.is_absolute() or not payload:
        raise PermissionError("Windows claim requires an absolute path and nonempty payload")
    advapi = ctypes.WinDLL("advapi32", use_last_error=True)
    kernel = ctypes.WinDLL("kernel32", use_last_error=True)
    sddl = f"D:P(A;;FA;;;SY)(A;;FA;;;BA)(A;;FA;;;{_user_sid()})"
    convert = advapi.ConvertStringSecurityDescriptorToSecurityDescriptorW
    convert.argtypes = [
        wintypes.LPCWSTR,
        wintypes.DWORD,
        ctypes.POINTER(ctypes.c_void_p),
        ctypes.c_void_p,
    ]
    convert.restype = wintypes.BOOL
    expected_sd = ctypes.c_void_p()
    if not convert(sddl, 1, ctypes.byref(expected_sd), None):
        raise ctypes.WinError(ctypes.get_last_error())
    handle = None
    held_identity: tuple[int, int] | None = None
    try:
        expected_dacl = _dacl_string(expected_sd)
        attributes = SecurityAttributes(ctypes.sizeof(SecurityAttributes), expected_sd, False)
        create = kernel.CreateFileW
        create.argtypes = [
            wintypes.LPCWSTR,
            wintypes.DWORD,
            wintypes.DWORD,
            ctypes.POINTER(SecurityAttributes),
            wintypes.DWORD,
            wintypes.DWORD,
            wintypes.HANDLE,
        ]
        create.restype = wintypes.HANDLE
        # CREATE_NEW is exclusive; write-through and FlushFileBuffers concern file
        # contents, without asserting crash durability of the directory entry.
        handle = create(str(path), 0xC0000000, 0, ctypes.byref(attributes), 1, 0x80000080, None)
        if handle == ctypes.c_void_p(-1).value:
            error = ctypes.get_last_error()
            if error in (80, 183):
                raise FileExistsError(error, "Attempt claim already exists", str(path))
            raise ctypes.WinError(error)
        written = wintypes.DWORD()
        write = kernel.WriteFile
        write.argtypes = [
            wintypes.HANDLE,
            ctypes.c_void_p,
            wintypes.DWORD,
            ctypes.POINTER(wintypes.DWORD),
            ctypes.c_void_p,
        ]
        write.restype = wintypes.BOOL
        buffer = ctypes.create_string_buffer(payload)
        if not write(
            handle, buffer, len(payload), ctypes.byref(written), None
        ) or written.value != len(payload):
            raise OSError("Incomplete Windows attempt-claim write")
        flush = kernel.FlushFileBuffers
        flush.argtypes = [wintypes.HANDLE]
        flush.restype = wintypes.BOOL
        if not flush(handle):
            raise ctypes.WinError(ctypes.get_last_error())
        position = kernel.SetFilePointerEx
        position.argtypes = [wintypes.HANDLE, ctypes.c_longlong, ctypes.c_void_p, wintypes.DWORD]
        position.restype = wintypes.BOOL
        if not position(handle, 0, None, 0):
            raise ctypes.WinError(ctypes.get_last_error())
        read = kernel.ReadFile
        read.argtypes = [
            wintypes.HANDLE,
            ctypes.c_void_p,
            wintypes.DWORD,
            ctypes.POINTER(wintypes.DWORD),
            ctypes.c_void_p,
        ]
        read.restype = wintypes.BOOL
        observed = ctypes.create_string_buffer(len(payload))
        count = wintypes.DWORD()
        if not read(handle, observed, len(payload), ctypes.byref(count), None):
            raise ctypes.WinError(ctypes.get_last_error())
        if count.value != len(payload) or observed.raw != payload:
            raise ValueError("Windows attempt-claim contents differ after flush")
        get_security = advapi.GetSecurityInfo
        get_security.argtypes = [
            wintypes.HANDLE,
            wintypes.DWORD,
            wintypes.DWORD,
            ctypes.c_void_p,
            ctypes.c_void_p,
            ctypes.c_void_p,
            ctypes.c_void_p,
            ctypes.POINTER(ctypes.c_void_p),
        ]
        get_security.restype = wintypes.DWORD
        actual_sd = ctypes.c_void_p()
        result = get_security(handle, 1, 4, None, None, None, None, ctypes.byref(actual_sd))
        if result:
            raise OSError(result, "Cannot read Windows attempt-claim ACL")
        try:
            if _dacl_string(actual_sd) != expected_dacl:
                raise PermissionError("Windows attempt-claim ACL differs from the approved policy")
        finally:
            kernel.LocalFree(actual_sd)
        held = path.stat()
        if not held.st_ino:
            raise ValueError("Windows attempt-claim file identity is unavailable")
        held_identity = (held.st_dev, held.st_ino)
    finally:
        if handle is not None and handle != ctypes.c_void_p(-1).value:
            kernel.CloseHandle.argtypes = [wintypes.HANDLE]
            kernel.CloseHandle(handle)
        kernel.LocalFree(expected_sd)
    observed_identity = path.stat()
    if held_identity != (observed_identity.st_dev, observed_identity.st_ino):
        raise ValueError("Windows attempt-claim path identity changed")
