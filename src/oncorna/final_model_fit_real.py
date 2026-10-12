"""Metadata-only real fitting preparation. No real fitting entry point is enabled."""

from __future__ import annotations

import hashlib
import json
import os
import subprocess
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml

from oncorna.final_fit_export import sequence_digest
from oncorna.final_fit_real import _create_local_claim
from oncorna.final_model_fit import (
    ATTEMPT_PATH,
    AUTHORIZATION_PATH,
    CONFIG_SHA256,
    FIT_EXPRESSION_SHA256,
    OUTPUT_DIRECTORY,
    SPLIT_SHA256,
    FitInputs,
    _canonical,
    _digest,
    expected_authorization,
    verify_authorization_record,
)
from oncorna.final_test import (
    FIT_COUNTS,
    FIT_EXPRESSION_PATH,
    AuthorizationError,
    RealSource,
    validate_protocol,
)
from oncorna.preprocessing import identity_digest, sha256_file


@dataclass(frozen=True)
class VerifiedFitBinding:
    """Values bound from frozen metadata and fitting-side true-label columns."""

    inputs: FitInputs
    config_sha256: str
    split_sha256: str
    expression_sha256: str
    test_count: int


@dataclass(frozen=True)
class VerifiedPreparation:
    binding: VerifiedFitBinding
    authorization_sha256: str
    authorization: dict[str, Any]
    revision: str


def _derive_binding(
    root: Path,
    config: dict[str, Any],
    *,
    synthetic_fixture: bool,
    class_counts: dict[str, int],
) -> VerifiedFitBinding:
    """Verify fitting artifact bytes and metadata without full-cohort or test access."""
    validate_protocol(config, synthetic_fixture=synthetic_fixture)
    source = RealSource(root)
    snapshot = source.preflight(config)
    train_count = config["partitions"]["original_train"]
    fit_count = config["partitions"]["final_fit_count"]
    test_count = config["partitions"]["final_test_count"]
    if (
        len(snapshot.fit_pairs) != fit_count
        or len(snapshot.test_pairs) != test_count
        or len(snapshot.gene_ids) != config["task"]["feature_count_before_filter"]
        or len(snapshot.fit_pairs[:train_count]) != train_count
        or len(snapshot.fit_pairs[train_count:]) != config["partitions"]["original_validation"]
        or len({patient for patient, _ in snapshot.fit_pairs + snapshot.test_pairs})
        != fit_count + test_count
        or len({sample for _, sample in snapshot.fit_pairs + snapshot.test_pairs})
        != fit_count + test_count
    ):
        raise ValueError("Frozen fitting, test or gene membership differs")
    labels = tuple(str(label) for label in source._fit_labels(snapshot))
    if len(labels) != fit_count or Counter(labels) != Counter(class_counts):
        raise ValueError("Fitting true-label counts differ from approved population")
    inputs = FitInputs(
        root / FIT_EXPRESSION_PATH,
        snapshot.fit_pairs,
        snapshot.gene_ids,
        labels,
        "true_label",
        {
            key: snapshot.input_hashes[key]
            for key in ("training_oof_predictions", "validation_predictions")
        },
    )
    return VerifiedFitBinding(
        inputs,
        sha256_file(root / "configs/final_test_v1.yaml"),
        snapshot.input_hashes["frozen_split"],
        snapshot.input_hashes["fit_expression"],
        test_count,
    )


def derive_real_fit_binding(root: Path) -> VerifiedFitBinding:
    """Bind the fixed 675/169 population without opening test-containing files."""
    root = root.resolve(strict=True)
    config_path = root / "configs/final_test_v1.yaml"
    if sha256_file(config_path) != CONFIG_SHA256:
        raise AuthorizationError("Frozen scientific configuration checksum differs")
    config = yaml.safe_load(config_path.read_bytes())
    binding = _derive_binding(root, config, synthetic_fixture=False, class_counts=FIT_COUNTS)
    if (
        binding.split_sha256 != SPLIT_SHA256
        or binding.expression_sha256 != FIT_EXPRESSION_SHA256
        or len(binding.inputs.fit_pairs) != 675
        or binding.test_count != 169
        or len(binding.inputs.gene_ids) != 20530
    ):
        raise AuthorizationError("Frozen fitting or expression identity differs")
    return binding


def _fixed_paths(root: Path, authorization_path: Path) -> tuple[Path, Path]:
    fixed = root / AUTHORIZATION_PATH
    marker = root / ATTEMPT_PATH
    output = root / OUTPUT_DIRECTORY
    if (
        os.path.normcase(os.path.abspath(authorization_path))
        != os.path.normcase(os.path.abspath(fixed))
        or fixed.parent.resolve(strict=True) != fixed.parent
        or fixed.resolve(strict=True) != fixed
        or marker.parent.resolve(strict=True) != marker.parent
        or output.parent.resolve(strict=True) != output.parent
    ):
        raise AuthorizationError("Fitting authorization or output path is redirected")
    if os.path.lexists(output) or os.path.lexists(marker):
        raise FileExistsError("A final-model fitting attempt already exists")
    return fixed, marker


def _revision(root: Path) -> str:
    top = subprocess.check_output(["git", "rev-parse", "--show-toplevel"], cwd=root, text=True)
    revision = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root, text=True).strip()
    dirty = subprocess.check_output(["git", "status", "--porcelain"], cwd=root, text=True)
    if _canonical(Path(top.strip())) != _canonical(root) or dirty.strip():
        raise AuthorizationError("Fitting requires the exact clean reviewed checkout")
    return revision


def verify_real_fit_preparation(root: Path, authorization_path: Path) -> VerifiedPreparation:
    """Compare local approval with independently derived, hash-verified inputs."""
    root = root.resolve(strict=True)
    fixed, _ = _fixed_paths(root, authorization_path)
    revision = _revision(root)
    binding = derive_real_fit_binding(root)
    raw = fixed.read_bytes()
    record = json.loads(raw)
    if not isinstance(record, dict):
        raise AuthorizationError("Fitting authorization must be an object")
    expected = expected_authorization(
        root,
        revision=revision,
        attempt_id=record.get("attempt_id"),
        approval_reference=record.get("approval_reference"),
        config_sha256=binding.config_sha256,
        split_sha256=binding.split_sha256,
        fit_expression_sha256=binding.expression_sha256,
        ordered_fit_identity_sha256=identity_digest(binding.inputs.fit_pairs),
        ordered_gene_list_sha256=sequence_digest(binding.inputs.gene_ids),
        ordered_true_label_sha256=_digest(binding.inputs.true_labels),
        label_table_hashes=binding.inputs.label_table_hashes,
    )
    verify_authorization_record(record, expected)
    return VerifiedPreparation(binding, hashlib.sha256(raw).hexdigest(), expected, revision)


def claim_real_fit_attempt(root: Path, authorization_path: Path) -> VerifiedPreparation:
    """Exclusively claim one prepared attempt; never fit or create model output."""
    root = root.resolve(strict=True)
    prepared = verify_real_fit_preparation(root, authorization_path)
    fixed, marker = _fixed_paths(root, authorization_path)
    if hashlib.sha256(fixed.read_bytes()).hexdigest() != prepared.authorization_sha256:
        raise AuthorizationError("Fitting authorization changed before claim")
    claim = {
        "attempt_id": prepared.authorization["attempt_id"],
        "approval_reference": prepared.authorization["approval_reference"],
        "authorization_sha256": prepared.authorization_sha256,
        "approved_commit": prepared.revision,
        "canonical_checkout": prepared.authorization["canonical_checkout"],
        "canonical_output": prepared.authorization["canonical_output"],
        "status": "claimed_before_fitting",
    }
    _create_local_claim(marker, claim)
    return prepared
