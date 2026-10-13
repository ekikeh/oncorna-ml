"""Metadata-only real fitting preparation. No real fitting entry point is enabled."""

from __future__ import annotations

import hashlib
import json
import os
import subprocess
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any, BinaryIO

import numpy as np
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
    _fit_and_seal,
    expected_authorization,
    verify_authorization_record,
    verify_sealed_fit,
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
    expression_stamp: tuple[int, int, int, int, int]


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
    expression = binding.inputs.expression_path
    if expression.resolve(strict=True) != expression:
        raise AuthorizationError("Fitting-only expression path is redirected")
    stamp = _stamp(expression.stat())
    if sha256_file(expression) != binding.expression_sha256 or _stamp(expression.stat()) != stamp:
        raise AuthorizationError("Fitting-only expression changed during preparation")
    return VerifiedPreparation(binding, hashlib.sha256(raw).hexdigest(), expected, revision, stamp)


def _claim_record(prepared: VerifiedPreparation) -> dict[str, str]:
    return {
        "attempt_id": prepared.authorization["attempt_id"],
        "approval_reference": prepared.authorization["approval_reference"],
        "authorization_sha256": prepared.authorization_sha256,
        "approved_commit": prepared.revision,
        "canonical_checkout": prepared.authorization["canonical_checkout"],
        "canonical_output": prepared.authorization["canonical_output"],
        "status": "claimed_before_fitting",
    }


def claim_real_fit_attempt(root: Path, authorization_path: Path) -> VerifiedPreparation:
    """Exclusively claim one prepared attempt; never fit or create model output."""
    root = root.resolve(strict=True)
    prepared = verify_real_fit_preparation(root, authorization_path)
    fixed, marker = _fixed_paths(root, authorization_path)
    if hashlib.sha256(fixed.read_bytes()).hexdigest() != prepared.authorization_sha256:
        raise AuthorizationError("Fitting authorization changed before claim")
    _create_local_claim(marker, _claim_record(prepared))
    return prepared


def _stamp(stat: os.stat_result) -> tuple[int, int, int, int, int]:
    """Available local file identity and high-resolution mutation metadata."""
    if not stat.st_dev or not stat.st_ino:
        raise PermissionError("Filesystem cannot provide a usable fitting-file identity")
    return (stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns)


def _hash_handle(handle: BinaryIO) -> str:
    handle.seek(0)
    digest = hashlib.sha256()
    for chunk in iter(lambda: handle.read(1024 * 1024), b""):
        digest.update(chunk)
    handle.seek(0)
    return digest.hexdigest()


def _verified_source_handle(
    handle: BinaryIO, path: Path, approved_sha256: str
) -> tuple[int, int, int, int, int]:
    if path.resolve(strict=True) != path:
        raise AuthorizationError("Fitting-only expression path is redirected")
    baseline = _stamp(os.fstat(handle.fileno()))
    if _stamp(path.stat()) != baseline or _hash_handle(handle) != approved_sha256:
        raise AuthorizationError("Fitting-only expression changed before fitting")
    if _stamp(os.fstat(handle.fileno())) != baseline or _stamp(path.stat()) != baseline:
        raise AuthorizationError("Fitting-only expression changed during verification")
    return baseline


def _run_claimed_fit(
    root: Path, authorization_path: Path | None, *, synthetic_fixture: bool
) -> dict[str, Any]:
    """Connect claim, verified inputs and common seal; real execution stays locked."""
    if not synthetic_fixture:
        raise AuthorizationError(
            "Real final-model fitting remains disabled pending separate approval"
        )
    root = root.resolve(strict=True)
    if synthetic_fixture and root == Path(__file__).resolve().parents[2]:
        raise PermissionError("Synthetic fitting cannot use the real project checkout")
    fixed = root / AUTHORIZATION_PATH
    if authorization_path is None:
        authorization_path = fixed
    prepared = verify_real_fit_preparation(root, authorization_path)
    inputs = prepared.binding.inputs
    if synthetic_fixture and (
        not inputs.fit_pairs
        or any(not item.startswith("SYN-") for pair in inputs.fit_pairs for item in pair)
        or inputs.gene_ids != tuple(f"G{i}" for i in range(6))
        or inputs.expression_path.resolve().is_relative_to(
            Path(__file__).resolve().parents[2] / "data"
        )
    ):
        raise PermissionError("Synthetic integration requires synthetic fitting-only inputs")
    claimed = claim_real_fit_attempt(root, authorization_path)
    if claimed != prepared:
        raise AuthorizationError("Fitting inputs or authorization changed before claim")
    marker = root / ATTEMPT_PATH
    marker_stamp = _stamp(marker.stat())
    claim_bytes = (
        json.dumps(_claim_record(claimed), sort_keys=True, separators=(",", ":")) + "\n"
    ).encode("utf-8")
    output = root / OUTPUT_DIRECTORY
    config = yaml.safe_load((root / "configs/final_test_v1.yaml").read_bytes())
    validate_protocol(config, synthetic_fixture=synthetic_fixture)
    expression = inputs.expression_path
    with expression.open("rb", buffering=0) as source:
        source_stamp = _verified_source_handle(
            source, expression, claimed.binding.expression_sha256
        )
        if source_stamp != claimed.expression_stamp:
            raise AuthorizationError("Fitting-only expression identity changed after preparation")

        def checkpoint(_stage: str) -> None:
            if (
                _revision(root) != claimed.revision
                or sha256_file(fixed) != claimed.authorization_sha256
                or marker.resolve(strict=True) != marker
                or _stamp(marker.stat()) != marker_stamp
                or marker.read_bytes() != claim_bytes
                or output.parent.resolve(strict=True) != output.parent
                or (output.exists() and output.resolve(strict=True) != output)
                or derive_real_fit_binding(root) != claimed.binding
                or _verified_source_handle(source, expression, claimed.binding.expression_sha256)
                != source_stamp
            ):
                raise AuthorizationError("Claimed fitting inputs, paths or revision changed")

        def validate_manifest(manifest: dict[str, Any]) -> None:
            expected = {
                "source_revision": claimed.revision,
                "config_sha256": claimed.authorization["config_sha256"],
                "split_sha256": claimed.authorization["split_sha256"],
                "fit_expression_sha256": claimed.authorization["fit_expression_sha256"],
                "fit_identity_sha256": claimed.authorization["ordered_fit_identity_sha256"],
                "original_gene_sha256": claimed.authorization["ordered_gene_list_sha256"],
                "ordered_true_label_sha256": claimed.authorization["ordered_true_label_sha256"],
                "label_table_hashes": claimed.authorization["label_table_hashes"],
                "class_order": claimed.authorization["class_order"],
                "model_parameters": claimed.authorization["model_parameters"],
                "attempt_id": claimed.authorization["attempt_id"],
                "authorization_reference": claimed.authorization["approval_reference"],
                "fit_count": len(inputs.fit_pairs),
                "original_gene_count": len(inputs.gene_ids),
            }
            if any(manifest.get(key) != value for key, value in expected.items()):
                raise AuthorizationError("Fitting manifest differs from the approved claim")

        checkpoint("before_values")
        values = np.load(source, allow_pickle=False)
        if (
            values.shape != (len(inputs.fit_pairs), len(inputs.gene_ids))
            or values.dtype != np.dtype("<f8")
            or not values.flags.c_contiguous
            or not np.isfinite(values).all()
        ):
            raise ValueError("Fitting-only expression layout or values differ")
        checkpoint("after_values")
        manifest = _fit_and_seal(
            inputs,
            config,
            output,
            values=values,
            artifact_sha256=claimed.binding.expression_sha256,
            config_sha256=claimed.binding.config_sha256,
            split_sha256=claimed.binding.split_sha256,
            source_revision=claimed.revision,
            attempt_id=claimed.authorization["attempt_id"],
            authorization_reference=claimed.authorization["approval_reference"],
            checkpoint=checkpoint,
            validate_manifest=validate_manifest,
        )
    validate_manifest(manifest)
    return verify_sealed_fit(output, manifest)
