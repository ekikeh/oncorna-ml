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


class FileIdInfo(ctypes.Structure):
    _fields_ = [
        ("VolumeSerialNumber", ctypes.c_ulonglong),
        ("FileId", ctypes.c_ubyte * 16),
    ]


def _file_identity(kernel: ctypes.WinDLL, handle: int) -> tuple[int, bytes]:
    """Query the volume and 128-bit file ID from this open Windows handle."""
    query = kernel.GetFileInformationByHandleEx
    query.argtypes = [wintypes.HANDLE, wintypes.DWORD, ctypes.c_void_p, wintypes.DWORD]
    query.restype = wintypes.BOOL
    info = FileIdInfo()
    if not query(handle, 0x12, ctypes.byref(info), ctypes.sizeof(info)):
        raise ctypes.WinError(ctypes.get_last_error())
    file_id = bytes(info.FileId)
    if not any(file_id):
        raise ValueError("Windows filesystem did not provide a usable file ID")
    return info.VolumeSerialNumber, file_id


def _close_handle(kernel: ctypes.WinDLL, handle: int) -> None:
    close = kernel.CloseHandle
    close.argtypes = [wintypes.HANDLE]
    close.restype = wintypes.BOOL
    if not close(handle):
        raise ctypes.WinError(ctypes.get_last_error())


def _free_descriptor(kernel: ctypes.WinDLL, descriptor: ctypes.c_void_p) -> None:
    free = kernel.LocalFree
    free.argtypes = [ctypes.c_void_p]
    free.restype = ctypes.c_void_p
    if free(descriptor):
        raise ctypes.WinError(ctypes.get_last_error())


def _path_identity(kernel: ctypes.WinDLL, path: Path) -> tuple[int, bytes]:
    """Open the fixed path afresh, without following a final reparse point."""
    open_file = kernel.CreateFileW
    open_file.argtypes = [
        wintypes.LPCWSTR,
        wintypes.DWORD,
        wintypes.DWORD,
        ctypes.c_void_p,
        wintypes.DWORD,
        wintypes.DWORD,
        wintypes.HANDLE,
    ]
    open_file.restype = wintypes.HANDLE
    handle = open_file(str(path), 0x80, 0, None, 3, 0x00200080, None)
    if handle == ctypes.c_void_p(-1).value:
        raise ctypes.WinError(ctypes.get_last_error())
    try:
        return _file_identity(kernel, handle)
    finally:
        _close_handle(kernel, handle)


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
        _free_descriptor(kernel, output)


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
    held_identity: tuple[int, bytes] | None = None
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
            _free_descriptor(kernel, actual_sd)
        held_identity = _file_identity(kernel, handle)
    finally:
        try:
            if handle is not None and handle != ctypes.c_void_p(-1).value:
                _close_handle(kernel, handle)
        finally:
            _free_descriptor(kernel, expected_sd)
    if held_identity != _path_identity(kernel, path):
        raise ValueError("Windows attempt-claim path identity changed")


def verify_existing_claim(path: Path, payload: bytes) -> None:
    """Recheck an existing claim's bytes, file identity and restrictive DACL."""
    if os.name != "nt" or not path.is_absolute() or not payload:
        raise PermissionError("Windows claim verification requires an absolute path")
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
    try:
        expected_dacl = _dacl_string(expected_sd)
        open_file = kernel.CreateFileW
        open_file.argtypes = [
            wintypes.LPCWSTR,
            wintypes.DWORD,
            wintypes.DWORD,
            ctypes.c_void_p,
            wintypes.DWORD,
            wintypes.DWORD,
            wintypes.HANDLE,
        ]
        open_file.restype = wintypes.HANDLE
        handle = open_file(str(path), 0x80020000, 0, None, 3, 0x00200080, None)
        if handle == ctypes.c_void_p(-1).value:
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
        observed = ctypes.create_string_buffer(len(payload) + 1)
        count = wintypes.DWORD()
        if not read(handle, observed, len(observed), ctypes.byref(count), None):
            raise ctypes.WinError(ctypes.get_last_error())
        if count.value != len(payload) or observed.raw[: count.value] != payload:
            raise ValueError("Windows attempt-claim contents differ")
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
                raise PermissionError("Windows attempt-claim ACL differs")
        finally:
            _free_descriptor(kernel, actual_sd)
        held_identity = _file_identity(kernel, handle)
    finally:
        try:
            if handle is not None and handle != ctypes.c_void_p(-1).value:
                _close_handle(kernel, handle)
        finally:
            _free_descriptor(kernel, expected_sd)
    if held_identity != _path_identity(kernel, path):
        raise ValueError("Windows attempt-claim path identity changed")
