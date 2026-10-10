"""Synthetic-only tests for exact-range final-fit export infrastructure."""

from __future__ import annotations

import copy
import json
import os
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

import oncorna.final_fit_export as export_module
from oncorna.final_fit_export import (
    ExportPlan,
    authorized_row_indices,
    export_synthetic,
    inspect_npy_header,
    plan_from_metadata,
    verify_completed_export,
    verify_manifest_contract,
    verify_source_header_at_test_access,
)
from oncorna.final_fit_real import (
    expected_authorization,
    run_real_export,
    verify_authorization_record,
)


def fixture_plan(tmp_path: Path, *, name: str = "result") -> tuple[ExportPlan, np.ndarray]:
    source = tmp_path / "source.npy"
    values = np.arange(32, dtype=np.float64).reshape(8, 4) + 0.25
    np.save(source, values)
    samples = [f"SYN-S-{i}" for i in range(8)]
    patients = [f"SYN-P-{i}" for i in range(8)]
    partitions = {
        "train": [7, 0, 4],
        "validation": [2, 5],
        "test": [1, 3, 6],
    }
    split = {
        "split_id": "split_v1",
        "split_unit": "patient_id",
        "source_cohort": {"path": "synthetic_cohort.tsv", "sha256": "c" * 64},
        "partitions": {
            name: {
                "n_patients": len(indices),
                "patient_ids": [patients[i] for i in indices],
                "sample_ids": [samples[i] for i in indices],
                "patient_sample_map": [
                    {"patient_id": patients[i], "sample_id": samples[i]} for i in indices
                ],
            }
            for name, indices in partitions.items()
        },
    }
    schema = {
        "sample_ids": samples,
        "gene_ids": [f"G{i}" for i in range(4)],
        "shape": [8, 4],
        "dtype": "float64",
        "orientation": "samples_by_genes",
        "expression_units": "log2(normalized_count + 1)",
    }
    split_path = tmp_path / "split.json"
    schema_path = tmp_path / "schema.json"
    split_path.write_text(json.dumps(split), encoding="utf-8")
    schema_path.write_text(json.dumps(schema), encoding="utf-8")
    plan = plan_from_metadata(
        source_array=source,
        split_path=split_path,
        schema_path=schema_path,
        output_dir=tmp_path / name,
        expected_split_sha256=export_module.sha256_file(split_path),
        expected_schema_sha256=export_module.sha256_file(schema_path),
        expected_all_gene_sha256="a" * 64,
        expected_source_sha256="b" * 64,
        expected_cohort_sha256="c" * 64,
        expected_fit_count=5,
        expected_test_count=3,
        expected_gene_count=4,
        approved_cohort_identifier="synthetic cohort",
        command="synthetic exact-range export",
    )
    return plan, values


def test_exact_order_boundary_rows_and_no_test_requests(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    plan, values = fixture_plan(tmp_path)
    requested: list[int] = []
    original = export_module._read_row

    def spy(handle: object, header: export_module.ArrayHeader, index: int) -> bytes:
        requested.append(index)
        return original(handle, header, index)

    monkeypatch.setattr(export_module, "_read_row", spy)
    manifest = export_synthetic(plan)
    selected = [7, 0, 4, 2, 5]
    assert requested == selected + selected
    assert set(requested).isdisjoint({1, 3, 6})
    actual = np.load(plan.output_dir / "fit_expression.npy", allow_pickle=False)
    assert actual.dtype == np.float64
    assert np.array_equal(actual, values[selected, :])
    assert manifest["shape"] == [5, 4]
    assert manifest["dtype_descr"] == "<f8"
    assert manifest["source_hash_verified_at_export"] is False
    header = inspect_npy_header(plan.source_array, (8, 4))
    assert verify_completed_export(plan.output_dir, plan, header) == manifest


@pytest.mark.parametrize("kind", ["fortran", "big_endian", "bad_magic", "truncated", "extra"])
def test_unsupported_or_corrupt_npy_is_rejected(tmp_path: Path, kind: str) -> None:
    plan, values = fixture_plan(tmp_path)
    if kind == "fortran":
        np.save(plan.source_array, np.asfortranarray(values))
    elif kind == "big_endian":
        np.save(plan.source_array, values.astype(">f8"))
    else:
        raw = plan.source_array.read_bytes()
        if kind == "bad_magic":
            raw = b"X" + raw[1:]
        elif kind == "truncated":
            raw = raw[:-8]
        else:
            raw += b"extra"
        plan.source_array.write_bytes(raw)
    with pytest.raises(ValueError):
        export_synthetic(plan)
    assert (plan.output_dir / "failure.json").exists()
    assert not (plan.output_dir / "success.json").exists()
    with pytest.raises(FileExistsError):
        export_synthetic(plan)


def test_membership_and_duplicate_indices_rejected_before_value_requests(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    plan, _ = fixture_plan(tmp_path)
    bad = copy.copy(plan)
    object.__setattr__(bad, "fit_pairs", plan.fit_pairs[:-1] + (plan.fit_pairs[0],))
    monkeypatch.setattr(export_module, "_read_row", lambda *_: pytest.fail("Value row read"))
    with pytest.raises(ValueError, match="duplicates"):
        authorized_row_indices(bad)
    with pytest.raises(ValueError):
        export_synthetic(bad)
    assert not (bad.output_dir / "success.json").exists()


def test_wrong_frozen_membership_rejected_by_metadata(tmp_path: Path) -> None:
    plan, _ = fixture_plan(tmp_path)
    split_path = tmp_path / "split.json"
    split = json.loads(split_path.read_text(encoding="utf-8"))
    split["partitions"]["train"]["patient_sample_map"][0]["sample_id"] = "SYN-UNKNOWN"
    split_path.write_text(json.dumps(split), encoding="utf-8")
    with pytest.raises(ValueError, match="approved population or schema"):
        plan_from_metadata(
            source_array=plan.source_array,
            split_path=split_path,
            schema_path=tmp_path / "schema.json",
            output_dir=tmp_path / "never",
            expected_split_sha256=export_module.sha256_file(split_path),
            expected_schema_sha256=export_module.sha256_file(tmp_path / "schema.json"),
            expected_all_gene_sha256="a" * 64,
            expected_source_sha256="b" * 64,
            expected_cohort_sha256="c" * 64,
            expected_fit_count=5,
            expected_test_count=3,
            expected_gene_count=4,
            approved_cohort_identifier="synthetic cohort",
            command="synthetic exact-range export",
        )


def test_manifest_version_fields_and_output_tampering(tmp_path: Path) -> None:
    plan, _ = fixture_plan(tmp_path)
    manifest = export_synthetic(plan)
    header = inspect_npy_header(plan.source_array, (8, 4))
    digest = export_module.sha256_file(plan.output_dir / "fit_expression.npy")
    for change in (
        {"version": "final_fit_v2"},
        {"unexpected": "value"},
        {"byte_order": "big"},
        {"repeat_selected_payload_sha256": "0" * 64},
    ):
        altered = {**manifest, **change}
        with pytest.raises(ValueError):
            verify_manifest_contract(altered, plan, header, digest)
    for field in export_module.MANIFEST_FIELDS:
        with pytest.raises(ValueError):
            verify_manifest_contract(
                {key: value for key, value in manifest.items() if key != field},
                plan,
                header,
                digest,
            )
    with (plan.output_dir / "fit_expression.npy").open("ab") as handle:
        handle.write(b"tampered")
    with pytest.raises(ValueError):
        verify_completed_export(plan.output_dir, plan, header)


def test_changed_fitting_source_values_change_artifact(tmp_path: Path) -> None:
    plan, values = fixture_plan(tmp_path)
    first = export_synthetic(plan)
    values[0, 0] += 10.0
    np.save(plan.source_array, values)
    second = export_synthetic(replace(plan, output_dir=tmp_path / "second"))
    assert first["fit_expression_sha256"] != second["fit_expression_sha256"]
    assert first["selected_payload_sha256"] != second["selected_payload_sha256"]


def test_identical_synthetic_exports_are_byte_deterministic(tmp_path: Path) -> None:
    plan, _ = fixture_plan(tmp_path)
    first = export_synthetic(plan)
    second_plan = replace(plan, output_dir=tmp_path / "second")
    second = export_synthetic(second_plan)
    assert first == second
    assert (plan.output_dir / "fit_expression.npy").read_bytes() == (
        second_plan.output_dir / "fit_expression.npy"
    ).read_bytes()


def test_interruption_preserves_attempt_and_refuses_retry(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    plan, _ = fixture_plan(tmp_path)
    original = export_module._read_row
    called = 0

    def interrupt(handle: object, header: export_module.ArrayHeader, index: int) -> bytes:
        nonlocal called
        called += 1
        if called == 2:
            raise KeyboardInterrupt("synthetic interruption")
        return original(handle, header, index)

    monkeypatch.setattr(export_module, "_read_row", interrupt)
    with pytest.raises(KeyboardInterrupt, match="synthetic interruption"):
        export_synthetic(plan)
    assert (plan.output_dir / "failure.json").exists()
    assert not (plan.output_dir / "success.json").exists()
    with pytest.raises(FileExistsError):
        export_synthetic(plan)


def test_interruption_before_success_does_not_seal(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    plan, _ = fixture_plan(tmp_path)
    original = export_module._atomic_json

    def interrupt(path: Path, value: dict[str, object]) -> None:
        if path.name == "success.json":
            raise KeyboardInterrupt("synthetic sealing interruption")
        original(path, value)

    monkeypatch.setattr(export_module, "_atomic_json", interrupt)
    with pytest.raises(KeyboardInterrupt, match="sealing interruption"):
        export_synthetic(plan)
    assert (plan.output_dir / "manifest.json").exists()
    assert not (plan.output_dir / "success.json").exists()
    assert (plan.output_dir / "failure.json").exists()
    with pytest.raises((FileNotFoundError, ValueError)):
        verify_completed_export(
            plan.output_dir, plan, inspect_npy_header(plan.source_array, (8, 4))
        )
    with pytest.raises(FileExistsError):
        export_synthetic(plan)


def test_real_id_plan_is_blocked_before_directory_creation(tmp_path: Path) -> None:
    plan, _ = fixture_plan(tmp_path)
    altered = replace(
        plan,
        fit_pairs=(("TCGA-XX-0001", plan.fit_pairs[0][1]),) + plan.fit_pairs[1:],
    )
    with pytest.raises(PermissionError, match="disabled"):
        export_synthetic(altered)
    assert not altered.output_dir.exists()


def test_source_modification_between_passes_preserves_failed_attempt(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    plan, _ = fixture_plan(tmp_path)
    original = export_module._read_row
    calls = 0

    def modify_after_first_pass(
        handle: object, header: export_module.ArrayHeader, index: int
    ) -> bytes:
        nonlocal calls
        row = original(handle, header, index)
        calls += 1
        if calls == len(plan.fit_pairs):
            modified = plan.source_array.stat().st_mtime_ns + 2_000_000_000
            os.utime(plan.source_array, ns=(modified, modified))
        return row

    monkeypatch.setattr(export_module, "_read_row", modify_after_first_pass)
    with pytest.raises(ValueError, match="metadata changed"):
        export_synthetic(plan)
    assert (plan.output_dir / "failure.json").exists()
    assert not (plan.output_dir / "success.json").exists()
    with pytest.raises(FileExistsError):
        export_synthetic(plan)


@pytest.mark.parametrize("field", ["inode", "size", "mtime_ns"])
def test_changed_source_path_identity_or_metadata_is_rejected(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, field: str
) -> None:
    plan, _ = fixture_plan(tmp_path)
    original = export_module._path_fingerprint
    calls = 0

    def changed_path(path: Path) -> export_module.FileFingerprint:
        nonlocal calls
        result = original(path)
        calls += 1
        if calls >= 3:
            return replace(result, **{field: getattr(result, field) + 1})
        return result

    monkeypatch.setattr(export_module, "_path_fingerprint", changed_path)
    with pytest.raises(ValueError, match="metadata changed"):
        export_synthetic(plan)
    assert (plan.output_dir / "failure.json").exists()
    assert not (plan.output_dir / "success.json").exists()


def test_strict_manifest_types_and_original_header_offset(tmp_path: Path) -> None:
    plan, _ = fixture_plan(tmp_path)
    manifest = export_synthetic(plan)
    header = inspect_npy_header(plan.source_array, (8, 4))
    digest = export_module.sha256_file(plan.output_dir / "fit_expression.npy")
    for altered in (
        {**manifest, "fit_patient_count": 5.0},
        {**manifest, "gene_count": True},
        {**manifest, "shape": [5.0, 4]},
        {**manifest, "source_array_shape": (8, 4)},
        {**manifest, "source_array_data_offset": True},
        {**manifest, "source_hash_verified_at_export": 0},
        {**manifest, "approved_cohort": {**manifest["approved_cohort"], "path": 1}},
        {**manifest, "selected_payload_sha256": 1},
        {**manifest, "software": {**manifest["software"], "numpy": 1}},
    ):
        with pytest.raises(ValueError):
            verify_manifest_contract(altered, plan, header, digest)
    verify_source_header_at_test_access(plan.source_array, (8, 4), header.data_offset)
    with pytest.raises(ValueError, match="offset"):
        verify_source_header_at_test_access(plan.source_array, (8, 4), header.data_offset + 16)


def test_missing_or_mismatched_frozen_provenance_and_count(tmp_path: Path) -> None:
    plan, _ = fixture_plan(tmp_path)
    original = plan_from_metadata
    arguments = {
        "source_array": plan.source_array,
        "split_path": tmp_path / "split.json",
        "schema_path": tmp_path / "schema.json",
        "output_dir": tmp_path / "never",
        "expected_split_sha256": plan.frozen_split_sha256,
        "expected_schema_sha256": plan.source_schema_sha256,
        "expected_all_gene_sha256": plan.source_all_gene_expression_sha256,
        "expected_source_sha256": plan.source_expression_sha256,
        "expected_cohort_sha256": plan.approved_cohort_sha256,
        "expected_fit_count": 5,
        "expected_test_count": 3,
        "expected_gene_count": 4,
        "approved_cohort_identifier": "synthetic cohort",
        "command": "synthetic exact-range export",
    }
    with pytest.raises(ValueError, match="checksum mismatch"):
        original(**{**arguments, "expected_split_sha256": "0" * 64})
    with pytest.raises(ValueError, match="population or schema"):
        original(**{**arguments, "expected_fit_count": 675})
    with pytest.raises(ValueError, match="SHA-256 provenance"):
        original(**{**arguments, "expected_source_sha256": ""})
    with pytest.raises(FileNotFoundError):
        original(**{**arguments, "schema_path": tmp_path / "missing-schema.json"})


def test_real_export_authorization_contract_remains_locked(tmp_path: Path) -> None:
    plan, _ = fixture_plan(tmp_path)
    record = expected_authorization(plan, "reviewed-synthetic-revision")
    verify_authorization_record(record, record)
    for altered in (
        {**record, "approved_commit": "different"},
        {**record, "manifest_version": "final_fit_v2"},
        {**record, "output_directory": "other"},
        {**record, "split_sha256": "0" * 64},
        {**record, "approved_commit": None},
        {**record, "unexpected": "value"},
    ):
        with pytest.raises(PermissionError):
            verify_authorization_record(altered, record)
    with pytest.raises(PermissionError, match="disabled"):
        run_real_export(tmp_path, tmp_path / "authorization.json")
    assert not (tmp_path / "data/processed/final_fit_v1").exists()
