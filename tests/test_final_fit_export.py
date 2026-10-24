"""Synthetic-only tests for exact-range final-fit export infrastructure."""

from __future__ import annotations

import copy
import json
import os
import uuid
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest
import yaml

import oncorna.final_fit_export as export_module
import oncorna.final_fit_real as real_module
import oncorna.final_test as final_test_module
import oncorna.windows_attempt_claim as claim_module
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


def test_interruption_after_sealing_is_not_reported_as_completed(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    plan, _ = fixture_plan(tmp_path)
    original = export_module._verify_export_seal

    def interrupt(*args: object, **kwargs: object) -> object:
        if kwargs.get("completed") is False:
            raise KeyboardInterrupt("before completion")
        return original(*args, **kwargs)

    monkeypatch.setattr(
        export_module,
        "_verify_export_seal",
        interrupt,
    )
    with pytest.raises(KeyboardInterrupt, match="before completion"):
        export_synthetic(plan)
    assert (plan.output_dir / "manifest.json").exists()
    assert (plan.output_dir / "success.json").exists()
    assert (plan.output_dir / "failure.json").exists()
    assert (
        '"state": "failed"'
        in (plan.output_dir / "export_state.jsonl").read_text(encoding="utf-8").splitlines()[-1]
    )
    with pytest.raises(FileExistsError):
        export_synthetic(plan)


def test_completed_append_failure_with_failed_append_failure_is_not_success(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    plan, _ = fixture_plan(tmp_path)
    original = export_module._append_state

    def fail_after_append(output: Path, state: str, **details: object) -> None:
        if state == "failed":
            raise OSError("synthetic failed-state write error")
        original(output, state, **details)
        if state == "completed":
            raise OSError("synthetic crash after completed append")

    monkeypatch.setattr(export_module, "_append_state", fail_after_append)
    with pytest.raises(OSError, match="after completed"):
        export_synthetic(plan)
    header = inspect_npy_header(plan.source_array, (8, 4))
    assert (plan.output_dir / "failure.json").exists()
    assert (
        '"state": "completed"'
        in (plan.output_dir / "export_state.jsonl").read_text(encoding="utf-8").splitlines()[-1]
    )
    with pytest.raises(ValueError, match="failed or contradictory"):
        verify_completed_export(plan.output_dir, plan, header)
    with pytest.raises(FileExistsError):
        export_synthetic(plan)


def test_completed_verifier_rejects_contradictory_history_and_tampered_seals(
    tmp_path: Path,
) -> None:
    plan, _ = fixture_plan(tmp_path)
    export_synthetic(plan)
    header = inspect_npy_header(plan.source_array, (8, 4))
    state = plan.output_dir / "export_state.jsonl"
    original_state = state.read_bytes()
    state.write_bytes(original_state + b'{"state":"failed"}\n')
    with pytest.raises(ValueError, match="failed or contradictory"):
        verify_completed_export(plan.output_dir, plan, header)
    state.write_bytes(original_state)
    manifest = plan.output_dir / "manifest.json"
    original_manifest = manifest.read_bytes()
    manifest.write_bytes(original_manifest + b" ")
    with pytest.raises(ValueError):
        verify_completed_export(plan.output_dir, plan, header)
    manifest.write_bytes(original_manifest)
    success = plan.output_dir / "success.json"
    success.write_bytes(success.read_bytes() + b" ")
    with pytest.raises(ValueError):
        verify_completed_export(plan.output_dir, plan, header)


def test_real_id_plan_is_blocked_before_directory_creation(tmp_path: Path) -> None:
    plan, _ = fixture_plan(tmp_path)
    altered = replace(
        plan,
        fit_pairs=(("TCGA-XX-0001", plan.fit_pairs[0][1]),) + plan.fit_pairs[1:],
    )
    with pytest.raises(PermissionError, match="synthetic"):
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


def synthetic_claim(marker: Path, record: dict[str, str]) -> None:
    """Test-only exclusive claim; never used by the real Windows entry point."""
    with marker.open("x", encoding="utf-8") as handle:
        json.dump(record, handle)


SYNTHETIC_ATTEMPT_ID = "61230c8b-bf74-4787-86a0-b8bf09ebc100"


def synthetic_authorization(plan: ExportPlan, revision: str, root: Path) -> dict[str, str]:
    return expected_authorization(
        plan,
        revision,
        root,
        attempt_id=SYNTHETIC_ATTEMPT_ID,
        approval_reference="SYNTHETIC-APPROVAL",
    )


def test_real_export_authorization_contract_remains_locked(tmp_path: Path) -> None:
    plan, _ = fixture_plan(tmp_path)
    plan = replace(plan, output_dir=tmp_path / real_module.OUTPUT_DIRECTORY)
    plan.output_dir.parent.mkdir(parents=True)
    record = synthetic_authorization(plan, "reviewed-synthetic-revision", tmp_path)
    verify_authorization_record(record, record)
    for altered in (
        {**record, "approved_commit": "different"},
        {**record, "attempt_id": str(uuid.uuid4())},
        {**record, "approval_reference": ""},
        {**record, "manifest_version": "final_fit_v2"},
        {**record, "output_directory": "other"},
        {**record, "split_sha256": "0" * 64},
        {**record, "approved_commit": None},
        {**record, "unexpected": "value"},
    ):
        with pytest.raises(PermissionError):
            verify_authorization_record(altered, record)
    with pytest.raises(PermissionError, match="UUIDv4"):
        expected_authorization(
            plan,
            "reviewed-synthetic-revision",
            tmp_path,
            attempt_id="not-a-uuid",
            approval_reference="SYNTHETIC-APPROVAL",
        )
    with pytest.raises(PermissionError, match="fixed local authorization path"):
        run_real_export(tmp_path, tmp_path / "authorization.json")
    assert not (tmp_path / "data/processed/final_fit_v1").exists()


def synthetic_real_preparation(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> tuple[Path, Path, dict[str, str]]:
    """Construct only synthetic metadata; the all-gene expression file stays absent."""
    root = tmp_path
    root.mkdir(parents=True, exist_ok=True)
    (root / "configs").mkdir()
    (root / "data/processed/preprocessing_v1").mkdir(parents=True)
    samples = [f"SYN-S-{i}" for i in range(844)]
    patients = [f"SYN-P-{i}" for i in range(844)]
    genes = [f"SYN-G-{i}" for i in range(20530)]
    sections = {"train": range(506), "validation": range(506, 675), "test": range(675, 844)}
    split = {
        "split_id": "split_v1",
        "split_unit": "patient_id",
        "source_cohort": {"path": "data/processed/classification_cohort.tsv", "sha256": "c" * 64},
        "partitions": {
            name: {
                "n_patients": len(indices),
                "patient_ids": [patients[index] for index in indices],
                "sample_ids": [samples[index] for index in indices],
                "patient_sample_map": [
                    {"patient_id": patients[index], "sample_id": samples[index]}
                    for index in indices
                ],
            }
            for name, indices in sections.items()
        },
    }
    schema = {
        "sample_ids": samples,
        "gene_ids": genes,
        "shape": [844, 20530],
        "dtype": "float64",
        "orientation": "samples_by_genes",
        "expression_units": "log2(normalized_count + 1)",
    }
    split_path = root / "data/processed/split_v1.json"
    schema_path = root / "data/processed/preprocessing_v1/all_gene_schema.json"
    split_path.write_text(json.dumps(split), encoding="utf-8")
    schema_path.write_text(json.dumps(schema), encoding="utf-8")
    split_hash = export_module.sha256_file(split_path)
    schema_hash = export_module.sha256_file(schema_path)
    config = {
        "inputs": {
            "all_gene_expression": "data/processed/preprocessing_v1/all_gene_expression.npy",
            "all_gene_schema": "data/processed/preprocessing_v1/all_gene_schema.json",
            "frozen_split": "data/processed/split_v1.json",
            "expected_sha256": {
                "all_gene_expression": "a" * 64,
                "all_gene_schema": schema_hash,
                "frozen_split": split_hash,
                "source_expression": "b" * 64,
                "cohort": "c" * 64,
            },
        }
    }
    config_path = root / "configs/final_test_v1.yaml"
    config_path.write_text(yaml.safe_dump(config), encoding="utf-8")
    monkeypatch.setattr(real_module, "CONFIG_SHA256", export_module.sha256_file(config_path))
    monkeypatch.setattr(real_module, "SPLIT_SHA256", split_hash)
    monkeypatch.setattr(final_test_module, "validate_protocol", lambda *_args: None)

    def git_read(args: list[str], **_kwargs: object) -> str:
        if args[-1] == "--show-toplevel":
            return str(root)
        return "synthetic-reviewed-revision\n" if args[-1] == "HEAD" else ""

    monkeypatch.setattr(real_module.subprocess, "check_output", git_read)
    auth_path = root / real_module.AUTHORIZATION_PATH
    plan = export_module.plan_from_metadata(
        source_array=root / config["inputs"]["all_gene_expression"],
        split_path=split_path,
        schema_path=schema_path,
        output_dir=root / real_module.OUTPUT_DIRECTORY,
        expected_split_sha256=split_hash,
        expected_schema_sha256=schema_hash,
        expected_all_gene_sha256="a" * 64,
        expected_source_sha256="b" * 64,
        expected_cohort_sha256="c" * 64,
        expected_fit_count=675,
        expected_test_count=169,
        expected_gene_count=20530,
        approved_cohort_identifier="TCGA-BRCA Xena primary-tumor cohort",
        command="one-attempt final_fit_v1 exact-range export",
    )
    record = synthetic_authorization(plan, "synthetic-reviewed-revision", root)
    auth_path.write_text(json.dumps(record), encoding="utf-8")
    return root, auth_path, record


def test_real_preparation_uses_only_frozen_synthetic_metadata(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    root, auth_path, _ = synthetic_real_preparation(tmp_path, monkeypatch)
    plan = real_module.prepare_real_export(root, auth_path)
    assert len(plan.fit_pairs) == 675
    assert len(plan.test_pairs) == 169
    assert len(plan.gene_ids) == 20530
    assert plan.output_dir == root / real_module.OUTPUT_DIRECTORY
    assert not plan.source_array.exists()
    assert not plan.output_dir.exists()


def test_real_preparation_rejects_revision_dirty_tree_or_wrong_record(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    root, auth_path, record = synthetic_real_preparation(tmp_path, monkeypatch)
    for change in (
        {"approved_commit": "wrong-revision"},
        {"authorization_schema": "final_fit_single_export_v1"},
        {"attempt_id": "not-a-uuid"},
        {"approval_reference": ""},
        {"output_directory": "data/processed/other"},
        {"ordered_fit_identity_sha256": "0" * 64},
        {"source_array_sha256": "0" * 64},
        {"manifest_version": "final_fit_v2"},
    ):
        auth_path.write_text(json.dumps({**record, **change}), encoding="utf-8")
        with pytest.raises(PermissionError):
            real_module.prepare_real_export(root, auth_path)
    auth_path.write_text(json.dumps(record), encoding="utf-8")

    def dirty_git(args: list[str], **_kwargs: object) -> str:
        if args[-1] == "--show-toplevel":
            return str(root)
        return "synthetic-reviewed-revision\n" if args[-1] == "HEAD" else " M file.py\n"

    monkeypatch.setattr(real_module.subprocess, "check_output", dirty_git)
    with pytest.raises(PermissionError, match="clean reviewed"):
        real_module.prepare_real_export(root, auth_path)


def test_real_preparation_rejects_wrong_frozen_membership(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    root, auth_path, _ = synthetic_real_preparation(tmp_path, monkeypatch)
    split_path = root / "data/processed/split_v1.json"
    split = json.loads(split_path.read_text(encoding="utf-8"))
    split["partitions"]["train"]["patient_sample_map"][0]["sample_id"] = "SYN-S-675"
    split_path.write_text(json.dumps(split), encoding="utf-8")
    monkeypatch.setattr(real_module, "SPLIT_SHA256", export_module.sha256_file(split_path))
    with pytest.raises(ValueError, match="population or schema"):
        real_module.prepare_real_export(root, auth_path)


def test_real_mode_denies_direct_shared_writer_before_output(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    plan, _ = fixture_plan(tmp_path)
    monkeypatch.setattr(export_module, "_read_row", lambda *_: pytest.fail("expression row read"))
    with pytest.raises(PermissionError, match="fixed authorization"):
        export_module._export_exact_ranges(plan, real=True)
    assert not plan.output_dir.exists()


def test_activated_path_extracts_only_synthetic_fit_rows_and_consumes_attempt(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    plan, values = fixture_plan(tmp_path)
    plan = replace(plan, output_dir=tmp_path / real_module.OUTPUT_DIRECTORY)
    plan.output_dir.parent.mkdir(parents=True)
    authorization_path = tmp_path / real_module.AUTHORIZATION_PATH
    authorization_path.write_text(
        json.dumps(synthetic_authorization(plan, "synthetic-reviewed-revision", tmp_path)),
        encoding="utf-8",
    )
    monkeypatch.setattr(real_module, "prepare_real_export", lambda *_args: plan)
    monkeypatch.setattr(
        real_module, "_checked_revision", lambda _root: "synthetic-reviewed-revision"
    )
    requested: list[int] = []
    original = export_module._read_row

    def observe(handle: object, header: export_module.ArrayHeader, index: int) -> bytes:
        requested.append(index)
        return original(handle, header, index)

    monkeypatch.setattr(export_module, "_read_row", observe)
    manifest = real_module.run_real_export(tmp_path, authorization_path)
    assert requested == [7, 0, 4, 2, 5] * 2
    assert set(requested).isdisjoint({1, 3, 6})
    assert np.array_equal(
        np.load(plan.output_dir / "fit_expression.npy", allow_pickle=False),
        values[[7, 0, 4, 2, 5], :],
    )
    assert manifest["source_hash_verified_at_export"] is False
    assert (plan.output_dir / "success.json").exists()
    claim = json.loads((tmp_path / real_module.ATTEMPT_MARKER_PATH).read_text(encoding="utf-8"))
    assert claim["policy"] == "final_fit_local_audit_v1"
    assert claim["attempt_id"] == SYNTHETIC_ATTEMPT_ID
    assert claim["authorization_sha256"] == export_module.sha256_file(authorization_path)
    assert claim["canonical_output"] == real_module._canonical_output(tmp_path)
    assert claim["config_sha256"] == real_module.CONFIG_SHA256
    assert claim["split_sha256"] == real_module.SPLIT_SHA256
    assert (
        verify_completed_export(
            plan.output_dir, plan, inspect_npy_header(plan.source_array, (8, 4))
        )
        == manifest
    )
    with pytest.raises(FileExistsError):
        real_module.run_real_export(tmp_path, authorization_path)


def test_activated_path_keeps_failed_attempt_and_refuses_retry(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    plan, _ = fixture_plan(tmp_path)
    plan = replace(plan, output_dir=tmp_path / real_module.OUTPUT_DIRECTORY)
    plan.output_dir.parent.mkdir(parents=True)
    authorization_path = tmp_path / real_module.AUTHORIZATION_PATH
    authorization_path.write_text(
        json.dumps(synthetic_authorization(plan, "synthetic-reviewed-revision", tmp_path)),
        encoding="utf-8",
    )
    monkeypatch.setattr(real_module, "prepare_real_export", lambda *_args: plan)
    monkeypatch.setattr(
        real_module, "_checked_revision", lambda _root: "synthetic-reviewed-revision"
    )
    original = export_module._read_row
    calls = 0

    def mutate(handle: object, header: export_module.ArrayHeader, index: int) -> bytes:
        nonlocal calls
        row = original(handle, header, index)
        calls += 1
        if calls == len(plan.fit_pairs):
            changed = plan.source_array.stat().st_mtime_ns + 2_000_000_000
            os.utime(plan.source_array, ns=(changed, changed))
        return row

    monkeypatch.setattr(export_module, "_read_row", mutate)
    with pytest.raises(ValueError, match="metadata changed"):
        real_module.run_real_export(tmp_path, authorization_path)
    assert (plan.output_dir / "failure.json").exists()
    assert not (plan.output_dir / "success.json").exists()
    assert (tmp_path / real_module.ATTEMPT_MARKER_PATH).exists()
    with pytest.raises(FileExistsError):
        real_module.run_real_export(tmp_path, authorization_path)


def test_export_record_cannot_authorize_final_test(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    root, authorization_path, _ = synthetic_real_preparation(tmp_path, monkeypatch)
    config_path = root / "configs/final_test_v1.yaml"
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    config["outputs"] = {"directory": "data/processed/final_test_v1"}
    with pytest.raises(final_test_module.AuthorizationError, match="does not bind"):
        final_test_module.verify_real_authorization(root, config_path, config, authorization_path)


def test_real_export_rejects_missing_record_and_existing_output(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    root, authorization_path, record = synthetic_real_preparation(tmp_path, monkeypatch)
    with pytest.raises(PermissionError, match="record is required"):
        real_module.run_real_export(root)
    authorization_path.unlink()
    with pytest.raises(FileNotFoundError):
        real_module.prepare_real_export(root, authorization_path)
    authorization_path.write_text(json.dumps(record), encoding="utf-8")
    output = root / real_module.OUTPUT_DIRECTORY
    output.mkdir(parents=True)
    with pytest.raises(FileExistsError, match="already exists"):
        real_module.prepare_real_export(root, authorization_path)


def test_authorization_cannot_be_copied_to_another_checkout(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    first, _, record = synthetic_real_preparation(tmp_path / "first", monkeypatch)
    second, second_auth, _ = synthetic_real_preparation(tmp_path / "second", monkeypatch)
    assert record["canonical_checkout"] != real_module._canonical_checkout(second)
    second_auth.write_text(json.dumps(record), encoding="utf-8")
    with pytest.raises(PermissionError, match="authorization"):
        real_module.prepare_real_export(second, second_auth)
    assert real_module._canonical_checkout(first / ".") == record["canonical_checkout"]
    assert real_module._canonical_output(first / ".") == record["canonical_output"]
    if os.name == "nt":
        assert (
            real_module._canonical_checkout(Path(str(first).swapcase()))
            == record["canonical_checkout"]
        )
        assert (
            real_module._canonical_output(Path(str(first).swapcase())) == record["canonical_output"]
        )


def test_redirected_output_parent_is_rejected(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    root, authorization_path, _ = synthetic_real_preparation(tmp_path, monkeypatch)
    original = root / "data/processed"
    alternate = root / "data/alternate_processed"
    original.rename(alternate)
    try:
        original.symlink_to(alternate, target_is_directory=True)
    except OSError:
        pytest.skip("Directory symlinks require a Windows privilege unavailable here")
    with pytest.raises(PermissionError):
        real_module.prepare_real_export(root, authorization_path)
    assert not (root / real_module.OUTPUT_DIRECTORY).exists()
    assert not (root / real_module.ATTEMPT_MARKER_PATH).exists()


def test_resolved_output_parent_mismatch_is_rejected_without_symlink_privilege(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    root, authorization_path, _ = synthetic_real_preparation(tmp_path, monkeypatch)
    fixed_parent = root / real_module.OUTPUT_DIRECTORY.parent
    original_resolve = Path.resolve

    def redirected(path: Path, *args: object, **kwargs: object) -> Path:
        if path == fixed_parent:
            return root / "data/elsewhere"
        return original_resolve(path, *args, **kwargs)

    monkeypatch.setattr(Path, "resolve", redirected)
    with pytest.raises(PermissionError, match="redirected"):
        real_module.prepare_real_export(root, authorization_path)


@pytest.mark.parametrize("failure", ["permission", "synchronization", "crash"])
def test_uncertain_claim_stays_consumed_and_blocks_source_access(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, failure: str
) -> None:
    plan, _ = fixture_plan(tmp_path)
    plan = replace(plan, output_dir=tmp_path / real_module.OUTPUT_DIRECTORY)
    plan.output_dir.parent.mkdir(parents=True)
    authorization_path = tmp_path / real_module.AUTHORIZATION_PATH
    authorization_path.write_text(
        json.dumps(synthetic_authorization(plan, "synthetic-reviewed-revision", tmp_path)),
        encoding="utf-8",
    )
    monkeypatch.setattr(real_module, "prepare_real_export", lambda *_args: plan)

    def git_read(args: list[str], **_kwargs: object) -> str:
        if args[-1] == "--show-toplevel":
            return str(tmp_path)
        return "synthetic-reviewed-revision\n" if args[-1] == "HEAD" else ""

    monkeypatch.setattr(real_module.subprocess, "check_output", git_read)
    monkeypatch.setattr(export_module, "_read_row", lambda *_: pytest.fail("source row accessed"))

    def uncertain(marker: Path, record: dict[str, str]) -> None:
        synthetic_claim(marker, record)
        raise OSError(f"simulated {failure} failure after exclusive claim")

    monkeypatch.setattr(real_module, "_create_local_claim", uncertain)
    with pytest.raises(OSError, match=failure):
        real_module.run_real_export(tmp_path, authorization_path)
    assert (tmp_path / real_module.ATTEMPT_MARKER_PATH).exists()
    assert not plan.output_dir.exists()
    with pytest.raises(FileExistsError):
        real_module.run_real_export(tmp_path, authorization_path)


def test_failure_before_local_claim_stops_before_source_access(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    plan, _ = fixture_plan(tmp_path)
    plan = replace(plan, output_dir=tmp_path / real_module.OUTPUT_DIRECTORY)
    plan.output_dir.parent.mkdir(parents=True)
    authorization_path = tmp_path / real_module.AUTHORIZATION_PATH
    authorization_path.write_text(
        json.dumps(synthetic_authorization(plan, "synthetic-reviewed-revision", tmp_path)),
        encoding="utf-8",
    )
    monkeypatch.setattr(real_module, "prepare_real_export", lambda *_args: plan)
    monkeypatch.setattr(
        real_module, "_checked_revision", lambda _root: "synthetic-reviewed-revision"
    )
    monkeypatch.setattr(
        real_module,
        "_create_local_claim",
        lambda *_: (_ for _ in ()).throw(KeyboardInterrupt("before claim")),
    )
    monkeypatch.setattr(export_module, "_read_row", lambda *_: pytest.fail("source row accessed"))
    with pytest.raises(KeyboardInterrupt, match="before claim"):
        real_module.run_real_export(tmp_path, authorization_path)
    assert not (tmp_path / real_module.ATTEMPT_MARKER_PATH).exists()
    assert not plan.output_dir.exists()


def test_windows_acl_failure_blocks_source_access(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    if os.name != "nt":
        pytest.skip("Windows-specific ACL creation path")
    plan, _ = fixture_plan(tmp_path)
    plan = replace(plan, output_dir=tmp_path / real_module.OUTPUT_DIRECTORY)
    plan.output_dir.parent.mkdir(parents=True)
    authorization_path = tmp_path / real_module.AUTHORIZATION_PATH
    authorization_path.write_text(
        json.dumps(synthetic_authorization(plan, "synthetic-reviewed-revision", tmp_path)),
        encoding="utf-8",
    )
    monkeypatch.setattr(real_module, "prepare_real_export", lambda *_args: plan)

    def git_read(args: list[str], **_kwargs: object) -> str:
        if args[-1] == "--show-toplevel":
            return str(tmp_path)
        return "synthetic-reviewed-revision\n" if args[-1] == "HEAD" else ""

    monkeypatch.setattr(real_module.subprocess, "check_output", git_read)
    monkeypatch.setattr(export_module, "_read_row", lambda *_: pytest.fail("source row accessed"))
    monkeypatch.setattr(
        claim_module,
        "create_and_verify_claim",
        lambda *_: (_ for _ in ()).throw(PermissionError("synthetic ACL failure")),
    )
    with pytest.raises(PermissionError, match="ACL failure"):
        real_module.run_real_export(tmp_path, authorization_path)
    assert not plan.output_dir.exists()


def test_direct_real_writer_rejects_a_caller_supplied_plan(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    plan, _ = fixture_plan(tmp_path)
    plan = replace(plan, output_dir=tmp_path / real_module.OUTPUT_DIRECTORY)
    plan.output_dir.parent.mkdir(parents=True)
    authorization_path = tmp_path / real_module.AUTHORIZATION_PATH
    authorization_path.write_text("{}", encoding="utf-8")
    monkeypatch.setattr(real_module, "prepare_real_export", lambda *_args: plan)
    monkeypatch.setattr(export_module, "_read_row", lambda *_: pytest.fail("source row accessed"))
    altered = replace(plan, fit_pairs=tuple(reversed(plan.fit_pairs)))
    with pytest.raises(PermissionError, match="differs from verified"):
        export_module._export_exact_ranges(
            altered, real=True, root=tmp_path, authorization_path=authorization_path
        )
    assert not (tmp_path / real_module.ATTEMPT_MARKER_PATH).exists()


def test_claim_rejects_modified_record_and_wrong_checkout_root(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    root, authorization_path, record = synthetic_real_preparation(tmp_path, monkeypatch)
    plan = real_module.prepare_real_export(root, authorization_path)
    authorization_path.write_text(
        json.dumps({**record, "canonical_output": "C:/wrong/output"}), encoding="utf-8"
    )
    with pytest.raises(PermissionError, match="authorization"):
        real_module.claim_real_export_attempt(root, authorization_path, plan)
    authorization_path.write_text(json.dumps(record), encoding="utf-8")

    def wrong_checkout(args: list[str], **_kwargs: object) -> str:
        if args[-1] == "--show-toplevel":
            return str(root.parent)
        return "synthetic-reviewed-revision\n" if args[-1] == "HEAD" else ""

    monkeypatch.setattr(real_module.subprocess, "check_output", wrong_checkout)
    with pytest.raises(PermissionError, match="checkout root"):
        real_module.claim_real_export_attempt(root, authorization_path, plan)
    assert not (root / real_module.ATTEMPT_MARKER_PATH).exists()
