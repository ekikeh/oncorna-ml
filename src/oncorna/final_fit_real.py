"""Metadata-only preparation for a future, separately authorized real export."""

from __future__ import annotations

import hashlib
import json
import os
import stat
import subprocess
import uuid
from datetime import datetime, timezone
from pathlib import Path

import yaml

from oncorna.final_fit_export import (
    VERSION,
    ExportPlan,
    _export_exact_ranges,
    plan_from_metadata,
    sequence_digest,
)
from oncorna.preprocessing import identity_digest, sha256_file

CONFIG_SHA256 = "1cadc4bce4d891f904a06409ba12794afa854abf8733276dc4484e2f074baa36"
SPLIT_SHA256 = "33e464aa58a69bed96740ef76b011ad23689acf2ed3baef79729b62057a2e26d"
OUTPUT_DIRECTORY = Path("data/processed/final_fit_v1")
AUTHORIZATION_PATH = Path("data/processed/final_fit_v1_authorization.json")
ATTEMPT_MARKER_PATH = Path("data/processed/final_fit_v1_attempt.json")


def _canonical_checkout(root: Path) -> str:
    return os.path.normcase(str(root.resolve(strict=True)))


def _canonical_output(root: Path) -> str:
    root = root.resolve(strict=True)
    fixed_parent = root / OUTPUT_DIRECTORY.parent
    if fixed_parent.resolve(strict=True) != fixed_parent:
        raise PermissionError("Real export output parent is redirected")
    return os.path.normcase(str(fixed_parent / OUTPUT_DIRECTORY.name))


def _require_fixed_authorization(root: Path, authorization_path: Path) -> Path:
    fixed = root / AUTHORIZATION_PATH
    if authorization_path.resolve() != fixed:
        raise PermissionError("Real export requires the fixed local authorization path")
    if fixed.parent.resolve(strict=True) != fixed.parent or fixed.resolve() != fixed:
        raise PermissionError("Real export authorization path is redirected")
    return fixed


def _checked_revision(root: Path) -> str:
    top = subprocess.check_output(["git", "rev-parse", "--show-toplevel"], cwd=root, text=True)
    if _canonical_checkout(Path(top.strip())) != _canonical_checkout(root):
        raise PermissionError("Real export requires the exact reviewed checkout root")
    revision = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root, text=True).strip()
    dirty = subprocess.check_output(["git", "status", "--porcelain"], cwd=root, text=True).strip()
    if dirty:
        raise PermissionError("Real export requires a clean reviewed Git revision")
    return revision


def expected_authorization(
    plan: ExportPlan,
    revision: str,
    root: Path,
    *,
    attempt_id: str,
    approval_reference: str,
) -> dict[str, str]:
    """Exact local-audit record; a human must approve each attempted extraction."""
    root = root.resolve(strict=True)
    if plan.output_dir != root / OUTPUT_DIRECTORY:
        raise PermissionError("Real export output differs from the fixed checkout path")
    try:
        parsed_id = uuid.UUID(attempt_id)
    except (TypeError, ValueError, AttributeError) as exc:
        raise PermissionError("Local export requires a canonical UUIDv4 attempt ID") from exc
    if parsed_id.version != 4 or str(parsed_id) != attempt_id:
        raise PermissionError("Local export requires a canonical UUIDv4 attempt ID")
    if (
        type(approval_reference) is not str
        or not 1 <= len(approval_reference) <= 200
        or not all(character.isprintable() for character in approval_reference)
    ):
        raise PermissionError("Local export requires a human approval reference")
    return {
        "authorization_schema": "final_fit_local_audit_v1",
        "intent": "approve_one_local_fitting_expression_attempt",
        "attempt_id": attempt_id,
        "approval_reference": approval_reference,
        "approved_commit": revision,
        "canonical_checkout": _canonical_checkout(root),
        "canonical_output": _canonical_output(root),
        "config_sha256": CONFIG_SHA256,
        "split_sha256": SPLIT_SHA256,
        "source_array_path": "data/processed/preprocessing_v1/all_gene_expression.npy",
        "source_array_sha256": plan.source_all_gene_expression_sha256,
        "schema_path": "data/processed/preprocessing_v1/all_gene_schema.json",
        "schema_sha256": plan.source_schema_sha256,
        "ordered_fit_identity_sha256": identity_digest(plan.fit_pairs),
        "ordered_gene_list_sha256": sequence_digest(plan.gene_ids),
        "output_directory": OUTPUT_DIRECTORY.as_posix(),
        "manifest_version": VERSION,
    }


def _create_local_claim(marker: Path, record: dict[str, str]) -> None:
    """Claim exclusively and verify contents; crash persistence remains unproven."""
    payload = (json.dumps(record, sort_keys=True, separators=(",", ":")) + "\n").encode("utf-8")
    if os.name == "nt":
        from oncorna.windows_attempt_claim import create_and_verify_claim

        create_and_verify_claim(marker, payload)
    else:
        descriptor = os.open(marker, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
        try:
            with os.fdopen(descriptor, "wb") as handle:
                os.fchmod(handle.fileno(), 0o600)
                if stat.S_IMODE(os.fstat(handle.fileno()).st_mode) != 0o600:
                    raise PermissionError("Attempt marker permissions are not restrictive")
                handle.write(payload)
                handle.flush()
                os.fsync(handle.fileno())
            directory = os.open(marker.parent, os.O_RDONLY | os.O_DIRECTORY)
            try:
                os.fsync(directory)
            finally:
                os.close(directory)
        except BaseException:
            # An incomplete claim remains visible and cannot be reused automatically.
            raise
    if marker.read_bytes() != payload:
        raise ValueError("Attempt claim contents changed after creation")


def _authorization_record(raw: bytes) -> dict[str, str]:
    record = json.loads(raw.decode("utf-8"))
    if not isinstance(record, dict):
        raise PermissionError("Local export authorization record is invalid")
    return record


def _expected_for_record(plan: ExportPlan, revision: str, root: Path, raw: bytes) -> dict[str, str]:
    record = _authorization_record(raw)
    return expected_authorization(
        plan,
        revision,
        root,
        attempt_id=record.get("attempt_id"),
        approval_reference=record.get("approval_reference"),
    )


def verify_authorization_record(record: object, expected: dict[str, str]) -> None:
    if (
        not isinstance(record, dict)
        or set(record) != set(expected)
        or any(
            type(record[key]) is not str or record[key] != value for key, value in expected.items()
        )
    ):
        raise PermissionError("Real fitting-only export authorization is absent or mismatched")


def prepare_real_export(root: Path, authorization_path: Path) -> ExportPlan:
    """Inspect frozen metadata only; never open the all-gene expression source."""
    from oncorna.final_test import validate_protocol

    root = root.resolve(strict=True)
    config_path = root / "configs/final_test_v1.yaml"
    split_path = root / "data/processed/split_v1.json"
    fixed_authorization = _require_fixed_authorization(root, authorization_path)
    output = root / OUTPUT_DIRECTORY
    _canonical_output(root)
    if output.exists():
        raise FileExistsError("A fitting-only export attempt already exists")
    if (root / ATTEMPT_MARKER_PATH).exists():
        raise FileExistsError("A one-time fitting-only export attempt was already claimed")
    config_bytes = config_path.read_bytes()
    if (
        hashlib.sha256(config_bytes).hexdigest() != CONFIG_SHA256
        or sha256_file(split_path) != SPLIT_SHA256
    ):
        raise ValueError("Frozen real-export configuration or split checksum mismatch")
    config = yaml.safe_load(config_bytes.decode("utf-8"))
    validate_protocol(config)
    expected = config["inputs"]["expected_sha256"]
    if (
        config["inputs"]["all_gene_expression"]
        != "data/processed/preprocessing_v1/all_gene_expression.npy"
        or config["inputs"]["all_gene_schema"]
        != "data/processed/preprocessing_v1/all_gene_schema.json"
        or config["inputs"]["frozen_split"] != "data/processed/split_v1.json"
    ):
        raise ValueError("Real export input paths differ from the frozen protocol")
    revision = _checked_revision(root)
    plan = plan_from_metadata(
        source_array=root / config["inputs"]["all_gene_expression"],
        split_path=split_path,
        schema_path=root / config["inputs"]["all_gene_schema"],
        output_dir=output,
        expected_split_sha256=SPLIT_SHA256,
        expected_schema_sha256=expected["all_gene_schema"],
        expected_all_gene_sha256=expected["all_gene_expression"],
        expected_source_sha256=expected["source_expression"],
        expected_cohort_sha256=expected["cohort"],
        expected_fit_count=675,
        expected_test_count=169,
        expected_gene_count=20530,
        approved_cohort_identifier="TCGA-BRCA Xena primary-tumor cohort",
        command="one-attempt final_fit_v1 exact-range export",
    )
    authorization_bytes = fixed_authorization.read_bytes()
    verify_authorization_record(
        _authorization_record(authorization_bytes),
        _expected_for_record(plan, revision, root, authorization_bytes),
    )
    return plan


def claim_real_export_attempt(root: Path, authorization_path: Path, plan: ExportPlan) -> None:
    """Consume the fixed one-attempt authorization before creating output files."""
    root = root.resolve(strict=True)
    if plan.output_dir != root / OUTPUT_DIRECTORY:
        raise PermissionError("Real export output path differs from the approved fixed path")
    _require_fixed_authorization(root, authorization_path)
    _canonical_output(root)
    revision = _checked_revision(root)
    authorization_bytes = authorization_path.read_bytes()
    verify_authorization_record(
        _authorization_record(authorization_bytes),
        _expected_for_record(plan, revision, root, authorization_bytes),
    )
    authorization = _authorization_record(authorization_bytes)
    marker = root / ATTEMPT_MARKER_PATH
    record = {
        "status": "attempt_claimed",
        "policy": "final_fit_local_audit_v1",
        "attempt_id": authorization["attempt_id"],
        "approval_reference": authorization["approval_reference"],
        "authorization_sha256": hashlib.sha256(authorization_bytes).hexdigest(),
        "approved_commit": revision,
        "canonical_checkout": _canonical_checkout(root),
        "canonical_output": _canonical_output(root),
        "config_sha256": CONFIG_SHA256,
        "split_sha256": SPLIT_SHA256,
        "source_array_sha256_recorded": plan.source_all_gene_expression_sha256,
        "schema_sha256": plan.source_schema_sha256,
        "ordered_fit_identity_sha256": identity_digest(plan.fit_pairs),
        "ordered_gene_list_sha256": sequence_digest(plan.gene_ids),
        "output_directory": OUTPUT_DIRECTORY.as_posix(),
        "claimed_at_utc": datetime.now(timezone.utc).isoformat(),
    }
    _create_local_claim(marker, record)


def run_real_export(root: Path, authorization_path: Path | None = None) -> dict[str, object]:
    """Prepare an approved local-audit attempt; never authorizes final-model execution."""
    if authorization_path is None:
        raise PermissionError("A separate real-export authorization record is required")
    plan = prepare_real_export(root, authorization_path)
    return _export_exact_ranges(plan, real=True, root=root, authorization_path=authorization_path)
