"""Sealed fitting-only stage. Real execution remains unconditionally disabled."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import platform
import re
import uuid
import warnings
from collections import Counter
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable

import joblib
import numpy as np
import pandas as pd
import sklearn
from sklearn.exceptions import ConvergenceWarning

from oncorna.final_fit_export import sequence_digest
from oncorna.final_test import LABELS, AuthorizationError, modeling_parameters, validate_protocol
from oncorna.modeling import make_pipeline
from oncorna.preprocessing import identity_digest, sha256_file

VERSION = "final_model_fit_v1"
AUTHORIZATION_VERSION = "final_model_fit_local_audit_v1"
CONFIG_SHA256 = "1cadc4bce4d891f904a06409ba12794afa854abf8733276dc4484e2f074baa36"
SPLIT_SHA256 = "33e464aa58a69bed96740ef76b011ad23689acf2ed3baef79729b62057a2e26d"
FIT_EXPRESSION_SHA256 = "54db5e3d8d931c22c19c5ac871a8da720f06c8f4f921759cb2d2031f3ae2382a"
OUTPUT_DIRECTORY = Path("data/processed/final_model_fit_v1")
AUTHORIZATION_PATH = Path("data/processed/final_model_fit_v1_authorization.json")
ATTEMPT_PATH = Path("data/processed/final_model_fit_v1_attempt.json")
STATES = ("reserved", "fitting", "fitted", "completed")
FILES = frozenset(
    {
        "pipeline.joblib",
        "training_gene_list.tsv",
        "training_manifest.json",
        "manifest.json",
        "success.json",
        "fit_state.jsonl",
    }
)
MANIFEST_FIELDS = frozenset(
    {
        "version",
        "scope",
        "source_revision",
        "config_sha256",
        "split_sha256",
        "fit_expression_sha256",
        "fit_identity_sha256",
        "original_gene_sha256",
        "ordered_true_label_sha256",
        "label_column",
        "label_table_hashes",
        "class_order",
        "model_parameters",
        "fit_count",
        "original_gene_count",
        "retained_gene_count",
        "retained_gene_sha256",
        "attempt_id",
        "authorization_reference",
        "authorization_sha256",
        "full_source_hash_verified",
        "source_row_equality_verified",
        "software",
        "output_hashes",
    }
)


@dataclass(frozen=True)
class FitInputs:
    """Ordered fitting-only values and ground-truth labels; no test-data fields."""

    expression_path: Path
    fit_pairs: tuple[tuple[str, str], ...]
    gene_ids: tuple[str, ...]
    true_labels: tuple[str, ...]
    label_column: str
    label_table_hashes: dict[str, str]


def _digest(values: tuple[str, ...]) -> str:
    return hashlib.sha256(
        json.dumps(values, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


def _sha(value: object) -> bool:
    return type(value) is str and re.fullmatch(r"[0-9a-f]{64}", value) is not None


def _exact_types(actual: Any, wanted: Any) -> bool:
    if type(actual) is not type(wanted):
        return False
    if isinstance(wanted, dict):
        return set(actual) == set(wanted) and all(
            _exact_types(actual[key], value) for key, value in wanted.items()
        )
    if isinstance(wanted, list):
        return len(actual) == len(wanted) and all(
            _exact_types(item, value) for item, value in zip(actual, wanted, strict=True)
        )
    return True


def _canonical(path: Path) -> str:
    return os.path.normcase(str(path.resolve(strict=True)))


def expected_authorization(
    root: Path,
    *,
    revision: str,
    attempt_id: str,
    approval_reference: str,
    config_sha256: str,
    split_sha256: str,
    fit_expression_sha256: str,
    ordered_fit_identity_sha256: str,
    ordered_gene_list_sha256: str,
    ordered_true_label_sha256: str,
    label_table_hashes: dict[str, str],
) -> dict[str, Any]:
    """Build the exact proposed record; this does not grant execution."""
    root = root.resolve(strict=True)
    parent = root / OUTPUT_DIRECTORY.parent
    if parent.resolve(strict=True) != parent:
        raise PermissionError("Fitting output parent is redirected")
    parsed = uuid.UUID(attempt_id)
    if parsed.version != 4 or str(parsed) != attempt_id:
        raise ValueError("A canonical UUIDv4 attempt ID is required")
    if re.fullmatch(r"[0-9a-f]{40}", revision) is None or any(
        not _sha(value)
        for value in (
            config_sha256,
            split_sha256,
            fit_expression_sha256,
            ordered_fit_identity_sha256,
            ordered_gene_list_sha256,
            ordered_true_label_sha256,
        )
    ):
        raise ValueError("Exact revision and SHA-256 bindings are required")
    if (
        type(approval_reference) is not str
        or not 1 <= len(approval_reference) <= 200
        or not all(character.isprintable() for character in approval_reference)
    ):
        raise ValueError("An identifiable human approval reference is required")
    if (
        type(label_table_hashes) is not dict
        or set(label_table_hashes) != {"training_oof_predictions", "validation_predictions"}
        or not all(_sha(value) for value in label_table_hashes.values())
    ):
        raise ValueError("Both fitting-label table hashes are required")
    return {
        "authorization_schema": AUTHORIZATION_VERSION,
        "intent": "approve_one_local_final_model_fit_attempt_only",
        "attempt_id": attempt_id,
        "approval_reference": approval_reference,
        "approved_commit": revision,
        "canonical_checkout": _canonical(root),
        "canonical_output": os.path.normcase(str(parent / OUTPUT_DIRECTORY.name)),
        "output_directory": OUTPUT_DIRECTORY.as_posix(),
        "config_sha256": config_sha256,
        "split_sha256": split_sha256,
        "fit_expression_sha256": fit_expression_sha256,
        "ordered_fit_identity_sha256": ordered_fit_identity_sha256,
        "ordered_gene_list_sha256": ordered_gene_list_sha256,
        "ordered_true_label_sha256": ordered_true_label_sha256,
        "label_table_hashes": dict(label_table_hashes),
        "class_order": list(LABELS),
        "model_parameters": {
            "expression_gt": 1.0,
            "minimum_sample_fraction": 0.20,
            "penalty": "l2",
            "C": 1.0,
            "solver": "lbfgs",
            "fit_intercept": True,
            "class_weight": None,
            "tol": 0.0001,
            "max_iter": 1000,
        },
        "manifest_version": VERSION,
    }


def verify_authorization_record(record: object, expected: dict[str, Any]) -> None:
    if not _exact_types(record, expected) or record != expected:
        raise AuthorizationError("Fitting authorization is missing or mismatched")


def verify_real_fit_authorization(root: Path, authorization_path: Path) -> str:
    """Verify a record against independently derived fitting metadata; never execute."""
    from oncorna.final_model_fit_real import verify_real_fit_preparation

    return verify_real_fit_preparation(root, authorization_path).authorization_sha256


def run_real_fit(root: Path, authorization_path: Path | None = None) -> dict[str, Any]:
    """F25 public entry point; the real branch remains unconditionally disabled."""
    from oncorna.final_model_fit_real import _run_claimed_fit

    return _run_claimed_fit(root, authorization_path, synthetic_fixture=False)


def _new_json(path: Path, value: dict[str, Any]) -> None:
    with path.open("x", encoding="utf-8", newline="\n") as handle:
        json.dump(value, handle, ensure_ascii=False, indent=2, sort_keys=True)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())


def _state(output: Path, state: str, **details: Any) -> None:
    with (output / "fit_state.jsonl").open("a", encoding="utf-8", newline="\n") as handle:
        handle.write(json.dumps({"state": state, **details}, sort_keys=True) + "\n")
        handle.flush()
        os.fsync(handle.fileno())


def _array_sha(value: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(value).tobytes()).hexdigest()


def _validate_inputs(inputs: FitInputs, config: dict[str, Any], artifact_sha256: str) -> np.ndarray:
    validate_protocol(config, synthetic_fixture=True)
    expected_count = config["partitions"]["final_fit_count"]
    if (
        inputs.label_column != "true_label"
        or len(inputs.fit_pairs) != expected_count
        or len(inputs.true_labels) != expected_count
        or len(inputs.gene_ids) != config["task"]["feature_count_before_filter"]
        or len(set(inputs.fit_pairs)) != expected_count
        or len({patient for patient, _ in inputs.fit_pairs}) != expected_count
        or len({sample for _, sample in inputs.fit_pairs}) != expected_count
        or len(set(inputs.gene_ids)) != len(inputs.gene_ids)
        or any(
            not identifier.startswith("SYN-") for pair in inputs.fit_pairs for identifier in pair
        )
        or any(not gene.startswith("SYN-") for gene in inputs.gene_ids)
        or set(inputs.true_labels) != set(LABELS)
        or set(inputs.label_table_hashes) != {"training_oof_predictions", "validation_predictions"}
        or any(not _sha(value) for value in inputs.label_table_hashes.values())
    ):
        raise ValueError("Synthetic fitting membership, schema or true-label provenance differs")
    if sha256_file(inputs.expression_path) != artifact_sha256:
        raise ValueError("Fitting-only expression artifact checksum differs")
    values = np.load(inputs.expression_path, mmap_mode="r", allow_pickle=False)
    if (
        values.shape != (expected_count, len(inputs.gene_ids))
        or values.dtype != np.dtype("<f8")
        or not values.flags.c_contiguous
        or not np.isfinite(values).all()
    ):
        raise ValueError("Fitting-only expression artifact differs from approved input")
    return values


def run_synthetic_fit(
    inputs: FitInputs,
    config: dict[str, Any],
    output: Path,
    *,
    artifact_sha256: str,
    config_sha256: str,
    split_sha256: str,
    expected_fit_identity_sha256: str,
    expected_gene_sha256: str,
    expected_true_label_sha256: str,
    source_revision: str = "synthetic-revision",
    attempt_id: str = "synthetic-attempt",
    authorization_reference: str = "synthetic-fixture",
) -> dict[str, Any]:
    """Exercise sealing with synthetic IDs and a synthetic array only."""
    project_data = Path(__file__).resolve().parents[2] / "data"
    if (
        inputs.expression_path.resolve().is_relative_to(project_data.resolve())
        or output.resolve().is_relative_to(project_data.resolve())
        or any(not _sha(value) for value in (artifact_sha256, config_sha256, split_sha256))
    ):
        raise PermissionError("Synthetic fitting cannot use project data or invalid hashes")
    if (
        identity_digest(inputs.fit_pairs) != expected_fit_identity_sha256
        or sequence_digest(inputs.gene_ids) != expected_gene_sha256
        or _digest(inputs.true_labels) != expected_true_label_sha256
    ):
        raise ValueError("Frozen fitting identity, gene or true-label order differs")
    values = _validate_inputs(inputs, config, artifact_sha256)
    return _fit_and_seal(
        inputs,
        config,
        output,
        values=values,
        artifact_sha256=artifact_sha256,
        config_sha256=config_sha256,
        split_sha256=split_sha256,
        source_revision=source_revision,
        attempt_id=attempt_id,
        authorization_reference=authorization_reference,
        execution_scope="synthetic",
        authorization_sha256=None,
        claimed_root=None,
        checkpoint=lambda _stage: None,
        validate_manifest=lambda _manifest: None,
    )


def _fit_and_seal(
    inputs: FitInputs,
    config: dict[str, Any],
    output: Path,
    *,
    values: np.ndarray,
    artifact_sha256: str,
    config_sha256: str,
    split_sha256: str,
    source_revision: str,
    attempt_id: str,
    authorization_reference: str,
    execution_scope: str,
    authorization_sha256: str | None,
    claimed_root: Path | None,
    checkpoint: Callable[[str], None],
    validate_manifest: Callable[[dict[str, Any]], None],
) -> dict[str, Any]:
    """Fit synthetic inputs or an independently checked synthetic claimed fixture."""
    if execution_scope == "synthetic":
        project_data = Path(__file__).resolve().parents[2] / "data"
        if (
            claimed_root is not None
            or authorization_sha256 is not None
            or output.resolve().is_relative_to(project_data.resolve())
            or inputs.expression_path.resolve().is_relative_to(project_data.resolve())
            or output.name == OUTPUT_DIRECTORY.name
            or not source_revision.startswith("synthetic-")
            or not attempt_id.startswith("synthetic-")
            or not authorization_reference.startswith("synthetic-")
        ):
            raise AuthorizationError("Synthetic fitting cannot claim production output or approval")
        _validate_inputs(inputs, config, artifact_sha256)
    elif execution_scope == "claimed_synthetic":
        if claimed_root is None or not _sha(authorization_sha256):
            raise AuthorizationError("A verified synthetic claim is required")
    else:
        raise AuthorizationError("Real fitting remains disabled pending separate approval")

    def checked(stage: str) -> None:
        if execution_scope == "claimed_synthetic":
            from oncorna.final_model_fit_real import verify_claimed_synthetic_context

            verify_claimed_synthetic_context(
                claimed_root,
                output,
                inputs,
                authorization_sha256,
                attempt_id,
                authorization_reference,
                source_revision,
                config_sha256,
                split_sha256,
                artifact_sha256,
            )
            approved = np.load(inputs.expression_path, allow_pickle=False)
            if (
                values.dtype != np.dtype("<f8")
                or not values.flags.c_contiguous
                or values.shape != approved.shape
                or not np.array_equal(values, approved)
            ):
                raise AuthorizationError("Fitting values differ from claimed expression artifact")
        checkpoint(stage)

    checked("before_output")
    if output.exists():
        raise FileExistsError("Fitting attempt already exists")
    output.mkdir(parents=True, exist_ok=False)
    stage = "reserved"
    try:
        _state(output, "reserved")
        _state(output, "fitting")
        stage = "fitting"
        checked("before_fit")
        frame = pd.DataFrame(
            np.asarray(values),
            index=tuple(sample for _, sample in inputs.fit_pairs),
            columns=inputs.gene_ids,
        )
        pipeline = make_pipeline(modeling_parameters(config), 1.0)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            pipeline.fit(frame, np.asarray(inputs.true_labels, dtype=object))
        checked("after_fit")
        if any(issubclass(item.category, ConvergenceWarning) for item in caught):
            raise RuntimeError("ConvergenceWarning: fitting cannot be sealed")
        gene_filter = pipeline.named_steps["gene_filter"]
        scaler = pipeline.named_steps["scaler"]
        classifier = pipeline.named_steps["classifier"]
        retained = tuple(gene_filter.selected_gene_ids_)
        required = math.ceil(0.20 * len(inputs.fit_pairs))
        if (
            not retained
            or gene_filter.fitting_sample_ids_ != tuple(frame.index)
            or gene_filter.required_count_ != required
            or gene_filter.n_samples_fit_ != len(inputs.fit_pairs)
            or np.any(np.asarray(scaler.n_samples_seen_) != len(inputs.fit_pairs))
            or tuple(classifier.classes_) != tuple(sorted(LABELS))
            or len(retained) != len(scaler.mean_)
        ):
            raise ValueError("Fitted preprocessing or classifier state differs from protocol")
        pending = output / ".pipeline.joblib.pending"
        with pending.open("xb") as handle:
            joblib.dump(pipeline, handle)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(pending, output / "pipeline.joblib")
        with (output / "training_gene_list.tsv").open("x", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["gene_id"])
            writer.writerows((gene,) for gene in retained)
            handle.flush()
            os.fsync(handle.fileno())
        training = {
            "fit_count": len(inputs.fit_pairs),
            "class_counts": dict(Counter(inputs.true_labels)),
            "class_order": list(LABELS),
            "prevalence_expression_gt": gene_filter.expression_gt,
            "prevalence_required_count": gene_filter.required_count_,
            "retained_gene_count": len(retained),
            "retained_gene_ids": list(retained),
            "retained_gene_sha256": sequence_digest(retained),
            "fitted_sample_count": gene_filter.n_samples_fit_,
            "scaler_sample_count": int(np.asarray(scaler.n_samples_seen_).min()),
            "scaler_mean_sha256": _array_sha(scaler.mean_),
            "scaler_scale_sha256": _array_sha(scaler.scale_),
            "classifier_classes": list(classifier.classes_),
            "classifier_n_iter": classifier.n_iter_.tolist(),
            "classifier_coef_sha256": _array_sha(classifier.coef_),
            "classifier_intercept_sha256": _array_sha(classifier.intercept_),
            "warnings": [
                {"category": item.category.__name__, "message": str(item.message)}
                for item in caught
            ],
        }
        _new_json(output / "training_manifest.json", training)
        _state(output, "fitted", retained_gene_count=len(retained))
        stage = "fitted"
        manifest = {
            "version": VERSION,
            "scope": execution_scope,
            "source_revision": source_revision,
            "config_sha256": config_sha256,
            "split_sha256": split_sha256,
            "fit_expression_sha256": artifact_sha256,
            "fit_identity_sha256": identity_digest(inputs.fit_pairs),
            "original_gene_sha256": sequence_digest(inputs.gene_ids),
            "ordered_true_label_sha256": _digest(inputs.true_labels),
            "label_column": "true_label",
            "label_table_hashes": dict(inputs.label_table_hashes),
            "class_order": list(LABELS),
            "model_parameters": {
                **modeling_parameters(config),
                "C": 1.0,
            },
            "fit_count": len(inputs.fit_pairs),
            "original_gene_count": len(inputs.gene_ids),
            "retained_gene_count": len(retained),
            "retained_gene_sha256": sequence_digest(retained),
            "attempt_id": attempt_id,
            "authorization_reference": authorization_reference,
            "authorization_sha256": authorization_sha256,
            "full_source_hash_verified": False,
            "source_row_equality_verified": False,
            "software": {
                "python": platform.python_version(),
                "numpy": np.__version__,
                "pandas": pd.__version__,
                "scikit_learn": sklearn.__version__,
                "joblib": joblib.__version__,
            },
            "output_hashes": {
                name: sha256_file(output / name)
                for name in ("pipeline.joblib", "training_gene_list.tsv", "training_manifest.json")
            },
        }
        validate_manifest(manifest)
        _new_json(output / "manifest.json", manifest)
        checked("before_seal")
        _new_json(
            output / "success.json",
            {"status": "sealed", "manifest_sha256": sha256_file(output / "manifest.json")},
        )
        verify_sealed_fit(output, manifest, completed=False)
        checked("before_complete")
        _state(
            output,
            "completed",
            manifest_sha256=sha256_file(output / "manifest.json"),
            success_sha256=sha256_file(output / "success.json"),
        )
        verify_sealed_fit(output, manifest)
        checked("after_complete")
        return manifest
    except BaseException as exc:
        try:
            _new_json(
                output / "failure.json",
                {
                    "stage": stage,
                    "error_type": type(exc).__name__,
                    "at_utc": datetime.now(timezone.utc).isoformat(),
                },
            )
            _state(output, "failed", stage=stage)
        except OSError:
            pass
        raise


def verify_sealed_fit(
    output: Path, expected: dict[str, Any] | None = None, *, completed: bool = True
) -> dict[str, Any]:
    """Verify seals and metadata without deserializing untrusted joblib bytes."""
    states = [
        json.loads(line)
        for line in (output / "fit_state.jsonl").read_text(encoding="utf-8").splitlines()
    ]
    wanted = list(STATES if completed else STATES[:-1])
    if (
        [row.get("state") for row in states] != wanted
        or (output / "failure.json").exists()
        or any(output.glob(".*.pending"))
        or {path.name for path in output.iterdir()} != FILES
    ):
        raise ValueError("Fitting run is incomplete or has contradictory state evidence")
    manifest = json.loads((output / "manifest.json").read_bytes())
    success = json.loads((output / "success.json").read_bytes())
    manifest_sha = sha256_file(output / "manifest.json")
    success_sha = sha256_file(output / "success.json")
    if success != {"status": "sealed", "manifest_sha256": manifest_sha}:
        raise ValueError("Fitting success marker does not seal manifest")
    if completed and (
        states[-1].get("manifest_sha256") != manifest_sha
        or states[-1].get("success_sha256") != success_sha
    ):
        raise ValueError("Fitting terminal state differs from sealed artifacts")
    if (
        manifest.get("version") != VERSION
        or set(manifest) != MANIFEST_FIELDS
        or manifest.get("scope") not in {"synthetic", "claimed_synthetic"}
        or manifest.get("class_order") != list(LABELS)
        or manifest.get("label_column") != "true_label"
        or manifest.get("full_source_hash_verified") is not False
        or manifest.get("source_row_equality_verified") is not False
        or type(manifest.get("output_hashes")) is not dict
        or set(manifest["output_hashes"])
        != {"pipeline.joblib", "training_gene_list.tsv", "training_manifest.json"}
    ):
        raise ValueError("Fitting manifest contract differs")
    if manifest["scope"] == "synthetic":
        if (
            manifest.get("authorization_sha256") is not None
            or output.name == OUTPUT_DIRECTORY.name
            or output.resolve().is_relative_to(
                (Path(__file__).resolve().parents[2] / "data").resolve()
            )
            or not all(
                isinstance(manifest.get(key), str) and manifest[key].startswith("synthetic-")
                for key in ("source_revision", "attempt_id", "authorization_reference")
            )
        ):
            raise ValueError("Synthetic seal cannot claim production authorization")
    else:
        if not _sha(manifest.get("authorization_sha256")):
            raise ValueError("Claimed seal lacks authorization checksum")
        from oncorna.final_model_fit_real import (
            derive_real_fit_binding,
            verify_claimed_synthetic_context,
        )

        root = output.parents[2]
        binding = derive_real_fit_binding(root)
        if (
            manifest["fit_identity_sha256"] != identity_digest(binding.inputs.fit_pairs)
            or manifest["original_gene_sha256"] != sequence_digest(binding.inputs.gene_ids)
            or manifest["ordered_true_label_sha256"] != _digest(binding.inputs.true_labels)
            or manifest["label_table_hashes"] != binding.inputs.label_table_hashes
            or manifest["fit_count"] != len(binding.inputs.fit_pairs)
            or manifest["original_gene_count"] != len(binding.inputs.gene_ids)
        ):
            raise ValueError("Claimed fitting manifest differs from frozen binding")
        verify_claimed_synthetic_context(
            root,
            output,
            binding.inputs,
            manifest["authorization_sha256"],
            manifest["attempt_id"],
            manifest["authorization_reference"],
            manifest["source_revision"],
            manifest["config_sha256"],
            manifest["split_sha256"],
            manifest["fit_expression_sha256"],
        )
    frozen_parameters = {
        "expression_gt": 1.0,
        "minimum_sample_fraction": 0.20,
        "penalty": "l2",
        "C": 1.0,
        "solver": "lbfgs",
        "fit_intercept": True,
        "class_weight": None,
        "tol": 0.0001,
        "max_iter": 1000,
    }
    if (
        not _exact_types(manifest.get("model_parameters"), frozen_parameters)
        or manifest["model_parameters"] != frozen_parameters
        or any(
            not _sha(manifest.get(key))
            for key in (
                "config_sha256",
                "split_sha256",
                "fit_expression_sha256",
                "fit_identity_sha256",
                "original_gene_sha256",
                "ordered_true_label_sha256",
                "retained_gene_sha256",
            )
        )
        or type(manifest.get("label_table_hashes")) is not dict
        or set(manifest["label_table_hashes"])
        != {"training_oof_predictions", "validation_predictions"}
        or any(not _sha(value) for value in manifest["label_table_hashes"].values())
        or type(manifest.get("fit_count")) is not int
        or type(manifest.get("original_gene_count")) is not int
        or type(manifest.get("retained_gene_count")) is not int
        or manifest["fit_count"] <= 0
        or manifest["original_gene_count"] <= 0
    ):
        raise ValueError("Fitting provenance or scientific parameters differ")
    if expected is not None and manifest != expected:
        raise ValueError("Fitting manifest differs from approved record")
    for name, digest in manifest["output_hashes"].items():
        if not _sha(digest) or sha256_file(output / name) != digest:
            raise ValueError(f"Fitting artifact checksum mismatch: {name}")
    training = json.loads((output / "training_manifest.json").read_bytes())
    with (output / "training_gene_list.tsv").open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    genes = [row[0] for row in rows[1:] if len(row) == 1]
    if (
        rows[0] != ["gene_id"]
        or len(genes) != len(rows) - 1
        or genes != training["retained_gene_ids"]
        or len(genes) != training["retained_gene_count"]
        or training["retained_gene_sha256"] != sequence_digest(tuple(genes))
        or manifest["retained_gene_sha256"] != training["retained_gene_sha256"]
        or manifest["retained_gene_count"] != len(genes)
        or manifest["fit_count"] != training["fit_count"]
        or training["class_order"] != list(LABELS)
        or training["prevalence_expression_gt"] != 1.0
        or training["prevalence_required_count"] != math.ceil(0.20 * manifest["fit_count"])
        or training["fitted_sample_count"] != manifest["fit_count"]
        or training["scaler_sample_count"] != manifest["fit_count"]
        or set(training["class_counts"]) != set(LABELS)
        or sum(training["class_counts"].values()) != manifest["fit_count"]
        or training["classifier_classes"] != sorted(LABELS)
    ):
        raise ValueError("Fitted gene list or training manifest differs")
    return manifest
