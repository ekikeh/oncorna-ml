"""Metadata-only preparation for a future, separately authorized real export."""

from __future__ import annotations

import hashlib
import json
import subprocess
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


def expected_authorization(plan: ExportPlan, revision: str) -> dict[str, str]:
    """Exact one-attempt record for the fixed output; approval must name this revision."""
    return {
        "authorization_schema": "final_fit_single_export_v1",
        "intent": "authorize_one_fitting_only_expression_export",
        "approved_commit": revision,
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

    root = root.resolve()
    config_path = root / "configs/final_test_v1.yaml"
    split_path = root / "data/processed/split_v1.json"
    fixed_authorization = root / AUTHORIZATION_PATH
    output = root / OUTPUT_DIRECTORY
    if authorization_path.resolve() != fixed_authorization.resolve():
        raise PermissionError("Real export requires the fixed local authorization path")
    if not output.resolve().is_relative_to(root):
        raise PermissionError("Real export output escapes the approved repository")
    if output.exists():
        raise FileExistsError("A fitting-only export attempt already exists")
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
    revision = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root, text=True).strip()
    dirty = subprocess.check_output(["git", "status", "--porcelain"], cwd=root, text=True).strip()
    if dirty:
        raise PermissionError("Real export requires a clean reviewed Git revision")
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
    record = json.loads(fixed_authorization.read_text(encoding="utf-8"))
    verify_authorization_record(record, expected_authorization(plan, revision))
    return plan


def run_real_export(root: Path, authorization_path: Path | None = None) -> dict[str, object]:
    """Prepare the fixed authorized plan; the shared writer still denies real mode."""
    if authorization_path is None:
        raise PermissionError("A separate real-export authorization record is required")
    plan = prepare_real_export(root, authorization_path)
    return _export_exact_ranges(plan, real=True)
