"""Exact-range, synthetic-only implementation of the final-fit export contract."""

from __future__ import annotations

import ast
import hashlib
import json
import os
import platform
import sys
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

from oncorna.preprocessing import identity_digest, sha256_file

VERSION = "final_fit_v1"
METHOD = "npy_c_float64_exact_row_ranges_v1"
MAGIC = b"\x93NUMPY"
MANIFEST_FIELDS = frozenset(
    {
        "version",
        "export_method",
        "source_expression_sha256",
        "source_all_gene_expression_sha256",
        "source_schema_sha256",
        "frozen_split_sha256",
        "approved_cohort",
        "ordered_fit_identity_sha256",
        "ordered_gene_list_sha256",
        "fit_patient_count",
        "gene_count",
        "shape",
        "dtype",
        "dtype_descr",
        "byte_order",
        "orientation",
        "expression_units",
        "source_array_shape",
        "source_array_data_offset",
        "source_array_layout",
        "source_hash_verified_at_export",
        "source_hash_verification_stage",
        "selected_payload_sha256",
        "repeat_selected_payload_sha256",
        "fit_expression_sha256",
        "software",
        "command",
    }
)


@dataclass(frozen=True)
class ArrayHeader:
    shape: tuple[int, int]
    data_offset: int
    row_bytes: int


@dataclass(frozen=True)
class FileFingerprint:
    device: int
    inode: int
    size: int
    mtime_ns: int
    ctime_ns: int


@dataclass(frozen=True)
class ExportPlan:
    source_array: Path
    output_dir: Path
    fit_pairs: tuple[tuple[str, str], ...]
    test_pairs: tuple[tuple[str, str], ...]
    sample_ids: tuple[str, ...]
    gene_ids: tuple[str, ...]
    frozen_split_sha256: str
    source_schema_sha256: str
    source_all_gene_expression_sha256: str
    source_expression_sha256: str
    approved_cohort_identifier: str
    approved_cohort_path: str
    approved_cohort_sha256: str
    expression_units: str
    command: str


def sequence_digest(values: tuple[str, ...]) -> str:
    data = json.dumps(values, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(data).hexdigest()


def _is_sha256(value: str) -> bool:
    return (
        isinstance(value, str)
        and len(value) == 64
        and all(character in "0123456789abcdef" for character in value)
    )


def _read_exact(handle: Any, count: int) -> bytes:
    pieces = []
    remaining = count
    while remaining:
        piece = handle.read(remaining)
        if not piece:
            raise ValueError("Truncated NumPy array")
        pieces.append(piece)
        remaining -= len(piece)
    return b"".join(pieces)


def _fingerprint(stat: os.stat_result) -> FileFingerprint:
    return FileFingerprint(
        stat.st_dev,
        stat.st_ino,
        stat.st_size,
        stat.st_mtime_ns,
        stat.st_ctime_ns,
    )


def _path_fingerprint(path: Path) -> FileFingerprint:
    return _fingerprint(path.stat())


def _assert_source_stable(handle: Any, path: Path, baseline: FileFingerprint) -> None:
    if _fingerprint(os.fstat(handle.fileno())) != baseline or _path_fingerprint(path) != baseline:
        raise ValueError("Source array identity, size, or modification metadata changed")


def inspect_npy_header_handle(handle: Any, expected_shape: tuple[int, int]) -> ArrayHeader:
    """Validate a NumPy header on the same descriptor used to extract rows."""
    if sys.byteorder != "little":
        raise ValueError("Exact-range exporter requires a little-endian host")
    handle.seek(0)
    if _read_exact(handle, 6) != MAGIC:
        raise ValueError("Invalid NumPy magic header")
    version = tuple(_read_exact(handle, 2))
    if version == (1, 0):
        header_length = int.from_bytes(_read_exact(handle, 2), "little")
    elif version == (2, 0):
        header_length = int.from_bytes(_read_exact(handle, 4), "little")
    else:
        raise ValueError("Unsupported NumPy array header version")
    if not 1 <= header_length <= 65536:
        raise ValueError("Invalid NumPy header length")
    try:
        header = ast.literal_eval(_read_exact(handle, header_length).decode("latin1").strip())
    except (SyntaxError, ValueError, UnicodeDecodeError) as exc:
        raise ValueError("Corrupt NumPy array header") from exc
    if not isinstance(header, dict) or set(header) != {"descr", "fortran_order", "shape"}:
        raise ValueError("Unexpected NumPy array header fields")
    if header["descr"] != "<f8" or header["fortran_order"] is not False:
        raise ValueError("Only C-order little-endian float64 arrays are supported")
    shape = header["shape"]
    if (
        not isinstance(shape, tuple)
        or len(shape) != 2
        or any(type(value) is not int or value <= 0 for value in shape)
        or shape != expected_shape
    ):
        raise ValueError("NumPy array shape differs from approved schema")
    offset = handle.tell()
    if offset % 16:
        raise ValueError("NumPy data offset is not 16-byte aligned")
    row_bytes = shape[1] * 8
    if os.fstat(handle.fileno()).st_size != offset + shape[0] * row_bytes:
        raise ValueError("NumPy file length differs from header and declared shape")
    return ArrayHeader(shape, offset, row_bytes)


def inspect_npy_header(path: Path, expected_shape: tuple[int, int]) -> ArrayHeader:
    """Read only the NumPy header and validate exact C-order float64 layout."""
    with path.open("rb", buffering=0) as handle:
        return inspect_npy_header_handle(handle, expected_shape)


def verify_source_header_at_test_access(
    path: Path, expected_shape: tuple[int, int], recorded_offset: int
) -> None:
    """Independently check the original header after full-source access is authorized."""
    if type(recorded_offset) is not int or recorded_offset < 0:
        raise ValueError("Invalid recorded source array data offset")
    if inspect_npy_header(path, expected_shape).data_offset != recorded_offset:
        raise ValueError("Original source array data offset differs from export manifest")


def plan_from_metadata(
    *,
    source_array: Path,
    split_path: Path,
    schema_path: Path,
    output_dir: Path,
    expected_split_sha256: str,
    expected_schema_sha256: str,
    expected_all_gene_sha256: str,
    expected_source_sha256: str,
    expected_cohort_sha256: str,
    expected_fit_count: int,
    expected_test_count: int,
    expected_gene_count: int,
    approved_cohort_identifier: str,
    command: str,
) -> ExportPlan:
    """Validate metadata and authorized rows before opening any expression value range."""
    if not all(
        _is_sha256(digest)
        for digest in (
            expected_split_sha256,
            expected_schema_sha256,
            expected_all_gene_sha256,
            expected_source_sha256,
            expected_cohort_sha256,
        )
    ):
        raise ValueError("Fitting export requires complete SHA-256 provenance")
    split_bytes = split_path.read_bytes()
    schema_bytes = schema_path.read_bytes()
    if hashlib.sha256(split_bytes).hexdigest() != expected_split_sha256:
        raise ValueError("Frozen split checksum mismatch")
    if hashlib.sha256(schema_bytes).hexdigest() != expected_schema_sha256:
        raise ValueError("All-gene schema checksum mismatch")
    split = json.loads(split_bytes.decode("utf-8"))
    schema = json.loads(schema_bytes.decode("utf-8"))
    if split.get("split_id") != "split_v1" or split.get("split_unit") != "patient_id":
        raise ValueError("Unexpected frozen split identity")
    parts = split["partitions"]
    fit = tuple(
        (item["patient_id"], item["sample_id"])
        for name in ("train", "validation")
        for item in parts[name]["patient_sample_map"]
    )
    test = tuple(
        (item["patient_id"], item["sample_id"]) for item in parts["test"]["patient_sample_map"]
    )
    samples = tuple(schema["sample_ids"])
    genes = tuple(schema["gene_ids"])
    if (
        len(fit) != expected_fit_count
        or len(test) != expected_test_count
        or len(genes) != expected_gene_count
        or len(samples) != len(fit) + len(test)
        or len({patient for patient, _ in fit + test}) != len(fit) + len(test)
        or len({sample for _, sample in fit + test}) != len(fit) + len(test)
        or len(set(samples)) != len(samples)
        or len(set(genes)) != len(genes)
        or set(samples) != {sample for _, sample in fit + test}
        or schema["shape"] != [len(samples), len(genes)]
        or schema["orientation"] != "samples_by_genes"
        or schema["dtype"] != "float64"
        or schema["expression_units"] != "log2(normalized_count + 1)"
        or split["source_cohort"]["sha256"] != expected_cohort_sha256
        or not approved_cohort_identifier
        or not command
    ):
        raise ValueError("Fitting export metadata differs from approved population or schema")
    for name in ("train", "validation", "test"):
        pairs = tuple(
            (item["patient_id"], item["sample_id"]) for item in parts[name]["patient_sample_map"]
        )
        if (
            len(pairs) != parts[name]["n_patients"]
            or tuple(parts[name]["patient_ids"]) != tuple(patient for patient, _ in pairs)
            or tuple(parts[name]["sample_ids"]) != tuple(sample for _, sample in pairs)
        ):
            raise ValueError(f"Frozen {name} membership fields disagree")
    return ExportPlan(
        source_array=source_array,
        output_dir=output_dir,
        fit_pairs=fit,
        test_pairs=test,
        sample_ids=samples,
        gene_ids=genes,
        frozen_split_sha256=expected_split_sha256,
        source_schema_sha256=expected_schema_sha256,
        source_all_gene_expression_sha256=expected_all_gene_sha256,
        source_expression_sha256=expected_source_sha256,
        approved_cohort_identifier=approved_cohort_identifier,
        approved_cohort_path=split["source_cohort"]["path"],
        approved_cohort_sha256=expected_cohort_sha256,
        expression_units=schema["expression_units"],
        command=command,
    )


def authorized_row_indices(plan: ExportPlan) -> tuple[int, ...]:
    positions = {sample: index for index, sample in enumerate(plan.sample_ids)}
    fit_samples = tuple(sample for _, sample in plan.fit_pairs)
    test_samples = {sample for _, sample in plan.test_pairs}
    if len(set(fit_samples)) != len(fit_samples) or set(fit_samples) & test_samples:
        raise ValueError("Fitting row request overlaps or duplicates a test row")
    try:
        indices = tuple(positions[sample] for sample in fit_samples)
    except KeyError as exc:
        raise ValueError("Fitting sample absent from approved schema") from exc
    if len(set(indices)) != len(indices):
        raise ValueError("Duplicate authorized row index")
    return indices


def _read_row(handle: Any, header: ArrayHeader, index: int) -> bytes:
    if not 0 <= index < header.shape[0]:
        raise ValueError("Requested row index outside NumPy array")
    handle.seek(header.data_offset + index * header.row_bytes)
    return _read_exact(handle, header.row_bytes)


def _atomic_json(path: Path, value: dict[str, Any]) -> None:
    temporary = path.with_name(f".{path.name}.pending")
    with temporary.open("x", encoding="utf-8", newline="\n") as handle:
        json.dump(value, handle, indent=2, sort_keys=True, ensure_ascii=False)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def _append_state(output: Path, state: str, **details: Any) -> None:
    path = output / "export_state.jsonl"
    with path.open("a", encoding="utf-8", newline="\n") as handle:
        handle.write(json.dumps({"state": state, **details}, sort_keys=True) + "\n")
        handle.flush()
        os.fsync(handle.fileno())


def expected_manifest(
    plan: ExportPlan, header: ArrayHeader, artifact_sha256: str
) -> dict[str, Any]:
    return {
        "version": VERSION,
        "export_method": METHOD,
        "source_expression_sha256": plan.source_expression_sha256,
        "source_all_gene_expression_sha256": plan.source_all_gene_expression_sha256,
        "source_schema_sha256": plan.source_schema_sha256,
        "frozen_split_sha256": plan.frozen_split_sha256,
        "approved_cohort": {
            "identifier": plan.approved_cohort_identifier,
            "path": plan.approved_cohort_path,
            "sha256": plan.approved_cohort_sha256,
        },
        "ordered_fit_identity_sha256": identity_digest(plan.fit_pairs),
        "ordered_gene_list_sha256": sequence_digest(plan.gene_ids),
        "fit_patient_count": len(plan.fit_pairs),
        "gene_count": len(plan.gene_ids),
        "shape": [len(plan.fit_pairs), len(plan.gene_ids)],
        "dtype": "float64",
        "dtype_descr": "<f8",
        "byte_order": "little",
        "orientation": "samples_by_genes",
        "expression_units": plan.expression_units,
        "source_array_shape": list(header.shape),
        "source_array_data_offset": header.data_offset,
        "source_array_layout": "C",
        "source_hash_verified_at_export": False,
        "source_hash_verification_stage": "authorized_test_access",
        "fit_expression_sha256": artifact_sha256,
        "software": {
            "exporter_version": VERSION,
            "python": platform.python_version(),
            "numpy": np.__version__,
        },
        "command": plan.command,
    }


def _exact_types(actual: Any, expected: Any) -> bool:
    if type(actual) is not type(expected):
        return False
    if isinstance(expected, dict):
        return set(actual) == set(expected) and all(
            _exact_types(actual[key], value) for key, value in expected.items()
        )
    if isinstance(expected, list):
        return len(actual) == len(expected) and all(
            _exact_types(item, value) for item, value in zip(actual, expected, strict=True)
        )
    return True


def verify_manifest_contract(
    manifest: dict[str, Any], plan: ExportPlan, header: ArrayHeader, artifact_sha256: str
) -> None:
    if manifest.get("version") != VERSION or set(manifest) != MANIFEST_FIELDS:
        raise ValueError("Unsupported or incomplete final-fit manifest version")
    expected = expected_manifest(plan, header, artifact_sha256)
    for key, value in expected.items():
        if key not in {"software", "command"}:
            if not _exact_types(manifest[key], value) or manifest[key] != value:
                raise ValueError(f"Final-fit manifest integrity mismatch: {key}")
    if (
        type(manifest["source_array_data_offset"]) is not int
        or manifest["source_array_data_offset"] < 0
        or manifest["source_array_data_offset"] % 16
    ):
        raise ValueError("Invalid original array data offset")
    for key in (
        "source_expression_sha256",
        "source_all_gene_expression_sha256",
        "source_schema_sha256",
        "frozen_split_sha256",
        "ordered_fit_identity_sha256",
        "ordered_gene_list_sha256",
        "fit_expression_sha256",
    ):
        if not _is_sha256(manifest[key]):
            raise ValueError(f"Invalid SHA-256 manifest field: {key}")
    if not _is_sha256(manifest["approved_cohort"]["sha256"]):
        raise ValueError("Invalid approved cohort SHA-256")
    software = manifest["software"]
    if (
        not isinstance(software, dict)
        or set(software) != {"exporter_version", "python", "numpy"}
        or software["exporter_version"] != VERSION
        or not isinstance(software["python"], str)
        or not software["python"].startswith("3.11.")
        or not isinstance(software["numpy"], str)
        or not software["numpy"]
        or not isinstance(manifest["command"], str)
        or not manifest["command"]
    ):
        raise ValueError("Final-fit export software or command provenance is incomplete")
    selected = manifest["selected_payload_sha256"]
    repeated = manifest["repeat_selected_payload_sha256"]
    if not _is_sha256(selected) or not _is_sha256(repeated) or selected != repeated:
        raise ValueError("Final-fit extraction repeat digest is missing or inconsistent")


def verify_completed_export(output: Path, plan: ExportPlan, header: ArrayHeader) -> dict[str, Any]:
    entries = [
        json.loads(line)
        for line in (output / "export_state.jsonl").read_text(encoding="utf-8").splitlines()
    ]
    if not entries or entries[-1]["state"] != "completed":
        raise ValueError("Fitting-only export has no terminal completed state")
    manifest_path = output / "manifest.json"
    success_path = output / "success.json"
    manifest_sha = sha256_file(manifest_path)
    if manifest_sha != entries[-1]["manifest_sha256"]:
        raise ValueError("Fitting-only manifest checksum differs from completed state")
    success_sha = sha256_file(success_path)
    if success_sha != entries[-1]["success_sha256"]:
        raise ValueError("Fitting-only success marker differs from completed state")
    if json.loads(success_path.read_text(encoding="utf-8")) != {
        "status": "sealed",
        "manifest_sha256": manifest_sha,
    }:
        raise ValueError("Fitting-only success marker does not bind manifest")
    artifact_sha = sha256_file(output / "fit_expression.npy")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    verify_manifest_contract(manifest, plan, header, artifact_sha)
    output_header = inspect_npy_header(
        output / "fit_expression.npy", (len(plan.fit_pairs), len(plan.gene_ids))
    )
    if output_header.row_bytes != header.row_bytes:
        raise ValueError("Fitting-only artifact row width differs from approved source")
    payload_digest = hashlib.sha256()
    with (output / "fit_expression.npy").open("rb", buffering=0) as handle:
        handle.seek(output_header.data_offset)
        for _ in plan.fit_pairs:
            payload_digest.update(_read_exact(handle, output_header.row_bytes))
    if payload_digest.hexdigest() != manifest["selected_payload_sha256"]:
        raise ValueError("Fitting-only artifact payload differs from selected source rows")
    return manifest


def _export_exact_ranges(plan: ExportPlan, *, real: bool) -> dict[str, Any]:
    """Shared row writer; a later reviewed code revision must unlock real mode."""
    if type(real) is not bool:
        raise PermissionError("Invalid fitting-only export mode")
    if real:
        raise PermissionError("Real fitting-only export remains disabled pending separate approval")
    else:
        all_ids = [value for pair in plan.fit_pairs + plan.test_pairs for value in pair]
        project_data = Path(__file__).resolve().parents[2] / "data"
        if (
            not all(value.startswith("SYN-") for value in all_ids)
            or plan.source_array.resolve().is_relative_to(project_data.resolve())
            or plan.output_dir.resolve().is_relative_to(project_data.resolve())
        ):
            raise PermissionError("Only synthetic fitting-only export is enabled")
    output = plan.output_dir
    output.mkdir(parents=True, exist_ok=False)
    stage = "reserved"
    try:
        _append_state(output, "reserved")
        indices = authorized_row_indices(plan)
        with plan.source_array.open("rb", buffering=0) as source:
            baseline = _fingerprint(os.fstat(source.fileno()))
            _assert_source_stable(source, plan.source_array, baseline)
            header = inspect_npy_header_handle(source, (len(plan.sample_ids), len(plan.gene_ids)))
            _assert_source_stable(source, plan.source_array, baseline)
            _append_state(output, "extracting", source_fingerprint=asdict(baseline))
            stage = "extracting"
            pending = output / ".fit_expression.npy.pending"
            selected_hash = hashlib.sha256()
            with pending.open("xb", buffering=0) as target:
                np.lib.format.write_array_header_1_0(
                    target,
                    {
                        "descr": "<f8",
                        "fortran_order": False,
                        "shape": (len(indices), len(plan.gene_ids)),
                    },
                )
                for index in indices:
                    row = _read_row(source, header, index)
                    selected_hash.update(row)
                    if target.write(row) != len(row):
                        raise OSError("Incomplete fitting-only row write")
                os.fsync(target.fileno())
            _assert_source_stable(source, plan.source_array, baseline)
            repeat_hash = hashlib.sha256()
            for index in indices:
                repeat_hash.update(_read_row(source, header, index))
            _assert_source_stable(source, plan.source_array, baseline)
            if selected_hash.digest() != repeat_hash.digest():
                raise ValueError("Repeated fitting-only extraction differs")
            _append_state(output, "source_stable", source_fingerprint=asdict(baseline))
        output_header = inspect_npy_header(pending, (len(indices), len(plan.gene_ids)))
        if output_header.row_bytes != header.row_bytes:
            raise ValueError("Exported NumPy row width differs from source")
        final_array = output / "fit_expression.npy"
        os.replace(pending, final_array)
        artifact_sha = sha256_file(final_array)
        manifest = expected_manifest(plan, header, artifact_sha)
        manifest["selected_payload_sha256"] = selected_hash.hexdigest()
        manifest["repeat_selected_payload_sha256"] = repeat_hash.hexdigest()
        verify_manifest_contract(manifest, plan, header, artifact_sha)
        _atomic_json(output / "manifest.json", manifest)
        manifest_sha = sha256_file(output / "manifest.json")
        _atomic_json(output / "success.json", {"status": "sealed", "manifest_sha256": manifest_sha})
        success_sha = sha256_file(output / "success.json")
        _append_state(output, "completed", manifest_sha256=manifest_sha, success_sha256=success_sha)
        verify_completed_export(output, plan, header)
        return manifest
    except BaseException as exc:
        try:
            _atomic_json(
                output / "failure.json",
                {
                    "stage": stage,
                    "error_type": type(exc).__name__,
                    "message": str(exc),
                    "at_utc": datetime.now(timezone.utc).isoformat(),
                },
            )
            _append_state(output, "failed", failed_stage=stage)
        except OSError:
            pass  # The exclusive directory remains an incomplete attempt.
        raise


def export_synthetic(plan: ExportPlan) -> dict[str, Any]:
    """Exercise the exact-range implementation only with synthetic fixture identities."""
    return _export_exact_ranges(plan, real=False)
