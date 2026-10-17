"""Single-attempt test-only consumer of the independently authenticated final fit."""

from __future__ import annotations

import csv
import hashlib
import json
import os
import re
import subprocess
import uuid
from pathlib import Path
from typing import Any, Callable

import joblib
import numpy as np
import yaml

from oncorna.final_fit_real import _create_local_claim
from oncorna.final_model_fit import CONFIG_SHA256 as FIT_CONFIG_SHA256
from oncorna.final_model_fit import verify_sealed_fit
from oncorna.final_model_fit_real import _verify_existing_claim
from oncorna.final_test import (
    LABELS,
    AuthorizationError,
    RealSource,
    RunState,
    Snapshot,
    Source,
    SyntheticFileSource,
    SyntheticFixtureSource,
    _atomic_write_json_new,
    _score,
    _verify_output_hashes,
    _verify_saved_predictions,
    _write_json_new,
    _write_predictions,
    validate_protocol,
    verify_completed_run,
)
from oncorna.preprocessing import identity_digest, sha256_file

# These identities are the operator-confirmed offline values, fixed in reviewed code.
# Real execution never accepts replacement values from an authorization file or CLI.
FIT_REVISION = "aa99b56653ccce0f83ae33bb5af2f454edfe9af5"
FIT_ATTEMPT = "6dd22770-3986-4d9e-a1f7-bc67f7041c26"
FIT_AUTH_SHA256 = "73b8a12c2fa8f1ca2c7ac4e2c331bad6bb9d555c99e499497a9ea9fc8345753e"
MODEL_SHA256 = "ea54f2254e72cbd865c737389a500406cb4ced52621c3f64c90176875a967f90"
MANIFEST_SHA256 = "151a40983477cdf63348fef4b70c1b2b874df9c49b5541017fa3087afbb88d53"
FIT_DIR = Path("data/processed/final_model_fit_v1")
OUTPUT_DIR = Path("data/processed/final_test_v1")
AUTH_PATH = Path("data/processed/final_test_v1_authorization.json")
CLAIM_PATH = Path("data/processed/final_test_v1_attempt.json")
CONFIG_PATH = Path("configs/final_test_v1.yaml")


def verify_approved_fit(root: Path) -> dict[str, Any]:
    """Authenticate bytes before any joblib load; inspect fitting evidence without loading it."""
    output = root / FIT_DIR
    if (
        sha256_file(output / "pipeline.joblib") != MODEL_SHA256
        or sha256_file(output / "manifest.json") != MANIFEST_SHA256
    ):
        raise AuthorizationError("Final model differs from operator-confirmed offline identities")
    manifest = verify_sealed_fit(output, allow_historical_revision=True)
    if (
        manifest["scope"] != "production"
        or manifest["source_revision"] != FIT_REVISION
        or manifest["attempt_id"] != FIT_ATTEMPT
        or manifest["authorization_sha256"] != FIT_AUTH_SHA256
        or manifest["output_hashes"]["pipeline.joblib"] != MODEL_SHA256
        or manifest["fit_count"] != 675
        or manifest["original_gene_count"] != 20530
        or manifest["retained_gene_count"] != 17623
    ):
        raise AuthorizationError("Sealed final fit differs from approved production attempt")
    return manifest


def expected_test_authorization(
    root: Path,
    revision: str,
    config_sha256: str,
    split_sha256: str,
    attempt_id: str,
    approval_reference: str,
) -> dict[str, str]:
    """Define the exact future approval record; constructing it grants no access."""
    parsed = uuid.UUID(attempt_id)
    if parsed.version != 4 or str(parsed) != attempt_id:
        raise ValueError("Test attempt ID must be a canonical UUIDv4")
    if re.fullmatch(r"[0-9a-f]{40}", revision) is None or any(
        re.fullmatch(r"[0-9a-f]{64}", value) is None for value in (config_sha256, split_sha256)
    ):
        raise ValueError("Exact source revision and checksums are required")
    if not approval_reference or not approval_reference.isprintable():
        raise ValueError("An identifiable human test approval is required")
    return {
        "authorization_schema": "final_test_only_single_attempt_v1",
        "intent": "authorize_one_test_only_evaluation",
        "attempt_id": attempt_id,
        "approval_reference": approval_reference,
        "approved_commit": revision,
        "canonical_checkout": os.path.normcase(str(root.resolve(strict=True))),
        "canonical_output": os.path.normcase(str(root / OUTPUT_DIR)),
        "output_directory": OUTPUT_DIR.as_posix(),
        "config_sha256": config_sha256,
        "split_sha256": split_sha256,
        "fit_revision": FIT_REVISION,
        "fit_attempt_id": FIT_ATTEMPT,
        "fit_authorization_sha256": FIT_AUTH_SHA256,
        "pipeline_sha256": MODEL_SHA256,
        "fit_manifest_sha256": MANIFEST_SHA256,
    }


def _claim_payload(record: dict[str, str], authorization_sha256: str) -> bytes:
    claim = {
        "attempt_id": record["attempt_id"],
        "approval_reference": record["approval_reference"],
        "approved_commit": record["approved_commit"],
        "authorization_sha256": authorization_sha256,
        "canonical_checkout": record["canonical_checkout"],
        "canonical_output": record["canonical_output"],
        "status": "claimed_before_test_access",
    }
    return (json.dumps(claim, sort_keys=True, separators=(",", ":")) + "\n").encode()


def verify_test_authorization(root: Path, authorization_path: Path) -> tuple[dict[str, str], str]:
    """Read a separately approved fixed-path record before claiming or test access."""
    root = root.resolve(strict=True)
    fixed = root / AUTH_PATH
    if authorization_path != fixed or fixed.resolve(strict=True) != fixed:
        raise AuthorizationError("Test authorization path is redirected or missing")
    if (root / OUTPUT_DIR).exists() or (root / CLAIM_PATH).exists():
        raise FileExistsError("Test attempt or output already exists")
    revision = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root, text=True).strip()
    if subprocess.check_output(["git", "status", "--porcelain"], cwd=root, text=True).strip():
        raise AuthorizationError("Test execution requires a clean reviewed checkout")
    config_sha = sha256_file(root / CONFIG_PATH)
    config = yaml.safe_load((root / CONFIG_PATH).read_bytes())
    validate_protocol(config)
    split_sha = sha256_file(root / config["inputs"]["frozen_split"])
    if split_sha != config["inputs"]["expected_sha256"]["frozen_split"]:
        raise AuthorizationError("Frozen split checksum differs")
    raw = fixed.read_bytes()
    record = json.loads(raw)
    if type(record) is not dict:
        raise AuthorizationError("Test authorization is not an object")
    expected = expected_test_authorization(
        root,
        revision,
        config_sha,
        split_sha,
        record.get("attempt_id"),
        record.get("approval_reference"),
    )
    if record != expected:
        raise AuthorizationError("Test approval does not bind this exact attempt and fit")
    return expected, hashlib.sha256(raw).hexdigest()


def verify_active_test_attempt(root: Path, required_state: str) -> dict[str, str]:
    """Recheck the production source boundary on every test-containing access."""
    root = root.resolve(strict=True)
    if root != Path(__file__).resolve().parents[2] or required_state not in {
        "test_access_started",
        "predictions_sealed",
    }:
        raise AuthorizationError("Production source access requires the reviewed checkout")
    fixed_paths = (root / AUTH_PATH, root / CLAIM_PATH, root / OUTPUT_DIR)
    try:
        canonical = all(path.resolve(strict=True) == path for path in fixed_paths)
    except OSError as exc:
        raise AuthorizationError("Test authorization, claim and output are required") from exc
    if not canonical:
        raise AuthorizationError("Test authorization, claim or output path is redirected")
    revision = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root, text=True).strip()
    if subprocess.check_output(["git", "status", "--porcelain"], cwd=root, text=True).strip():
        raise AuthorizationError("Test source access requires a clean reviewed checkout")
    config_path = root / CONFIG_PATH
    config_sha = sha256_file(config_path)
    if config_sha != FIT_CONFIG_SHA256:
        raise AuthorizationError("Frozen scientific configuration differs")
    config = yaml.safe_load(config_path.read_bytes())
    validate_protocol(config)
    for key, fixed in (
        ("source_expression", "data/raw/TCGA-BRCA_HiSeqV2.tsv.gz"),
        ("all_gene_expression", "data/processed/preprocessing_v1/all_gene_expression.npy"),
        ("cohort", "data/processed/classification_cohort.tsv"),
    ):
        path = root / fixed
        if config["inputs"][key] != fixed or path.resolve(strict=True) != path:
            raise AuthorizationError(f"Production source path differs: {key}")
    split_sha = sha256_file(root / config["inputs"]["frozen_split"])
    if split_sha != config["inputs"]["expected_sha256"]["frozen_split"]:
        raise AuthorizationError("Frozen split differs before source access")
    raw = (root / AUTH_PATH).read_bytes()
    record = json.loads(raw)
    if type(record) is not dict:
        raise AuthorizationError("Test approval record differs")
    expected = expected_test_authorization(
        root,
        revision,
        config_sha,
        split_sha,
        record.get("attempt_id"),
        record.get("approval_reference"),
    )
    if record != expected:
        raise AuthorizationError("Test approval differs before source access")
    authorization_sha256 = hashlib.sha256(raw).hexdigest()
    _verify_existing_claim(root / CLAIM_PATH, _claim_payload(expected, authorization_sha256))
    entries = [
        json.loads(line)
        for line in (root / OUTPUT_DIR / "run_state.jsonl").read_text(encoding="utf-8").splitlines()
    ]
    wanted = ["reserved", "test_access_started"]
    if required_state == "predictions_sealed":
        wanted.append("predictions_sealed")
    if (
        [entry.get("state") for entry in entries] != wanted
        or [entry.get("sequence") for entry in entries] != list(range(1, len(wanted) + 1))
        or any(
            entry.get("attempt_id") != expected["attempt_id"]
            or entry.get("authorization_sha256") != authorization_sha256
            or entry.get("source_revision") != revision
            or entry.get("model_sha256") != MODEL_SHA256
            for entry in entries
        )
    ):
        raise AuthorizationError("Active test audit state differs from claimed attempt")
    return {**expected, "authorization_sha256": authorization_sha256}


def _validate_pipeline(pipeline: Any, snapshot: Snapshot, fit_manifest: dict[str, Any]) -> Any:
    if list(pipeline.named_steps) != ["gene_filter", "scaler", "classifier"]:
        raise ValueError("Loaded pipeline has unexpected steps")
    gene_filter = pipeline.named_steps["gene_filter"]
    scaler = pipeline.named_steps["scaler"]
    classifier = pipeline.named_steps["classifier"]
    if (
        tuple(gene_filter.feature_names_in_) != snapshot.gene_ids
        or tuple(gene_filter.selected_gene_ids_) != tuple(fit_manifest["retained_gene_ids"])
        or gene_filter.n_samples_fit_ != fit_manifest["fit_count"]
        or gene_filter.required_count_ != fit_manifest["prevalence_required_count"]
        or np.any(np.asarray(scaler.n_samples_seen_) != fit_manifest["fit_count"])
        or tuple(classifier.classes_) != tuple(fit_manifest["classifier_classes"])
        or set(classifier.classes_) != set(LABELS)
    ):
        raise ValueError("Loaded pipeline differs from sealed fitting state")
    return classifier


def verify_completed_test_only(output: Path) -> dict[str, Any]:
    """Independently reconcile the prediction seal with completed run evidence."""
    manifest = verify_completed_run(output)
    entries = [json.loads(line) for line in (output / "run_state.jsonl").read_text().splitlines()]
    if (
        [entry["state"] for entry in entries]
        != ["reserved", "test_access_started", "predictions_sealed", "completed"]
        or (output / "failure.json").exists()
        or any(output.glob(".*.pending"))
    ):
        raise ValueError("Test-only state history is incomplete or contradictory")
    seal = json.loads((output / "prediction_seal.json").read_bytes())
    prediction_sha = sha256_file(output / "test_predictions.tsv")
    identity_sha = manifest["test_identity_sha256"]
    binding = {
        "attempt_id": manifest["attempt_id"],
        "authorization_sha256": manifest["authorization_sha256"],
        "source_revision": manifest["source_revision"],
        "model_sha256": manifest["model_sha256"],
    }
    if (
        seal
        != {
            "status": "sealed_before_labels",
            **binding,
            "test_identity_sha256": identity_sha,
            "prediction_sha256": prediction_sha,
            "row_count": manifest["test_count"],
            "class_order": list(LABELS),
        }
        or any(any(entry.get(key) != value for key, value in binding.items()) for entry in entries)
        or entries[2]["prediction_sha256"] != prediction_sha
        or entries[2].get("test_identity_sha256") != identity_sha
        or entries[2]["seal_sha256"] != sha256_file(output / "prediction_seal.json")
        or manifest["sealed_prediction_sha256"] != prediction_sha
        or manifest["output_hashes"]["test_predictions.tsv"] != prediction_sha
    ):
        raise ValueError("Prediction seal, audit state or final manifest differs")
    with (output / "test_predictions.tsv").open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = [
            "patient_id",
            "sample_id",
            "predicted_label",
            *(f"probability_{label}" for label in LABELS),
        ]
        if reader.fieldnames != fields:
            raise ValueError("Sealed prediction schema differs")
        rows = list(reader)
    if (
        len(rows) != manifest["test_count"]
        or len({row["patient_id"] for row in rows}) != len(rows)
        or len({row["sample_id"] for row in rows}) != len(rows)
        or identity_digest(tuple((row["patient_id"], row["sample_id"]) for row in rows))
        != identity_sha
    ):
        raise ValueError("Sealed prediction identities differ")
    for row in rows:
        probabilities = np.asarray([float(row[f"probability_{label}"]) for label in LABELS])
        if (
            row["predicted_label"] not in LABELS
            or not np.isfinite(probabilities).all()
            or np.any(probabilities < 0)
            or np.any(probabilities > 1)
            or not np.isclose(probabilities.sum(), 1, atol=1e-9, rtol=0)
        ):
            raise ValueError("Sealed prediction probabilities differ")
    return manifest


def _verify_claimed_real_context(
    source: RealSource,
    output: Path,
    pipeline_path: Path,
    fit_manifest: dict[str, Any],
    config: dict[str, Any],
    authorization_sha256: str,
    model_sha256: str,
    source_revision: str,
) -> dict[str, str]:
    """Guard the shared executor even when it is called outside the CLI."""
    root = source.root.resolve(strict=True)
    if (
        output != root / OUTPUT_DIR
        or pipeline_path != root / FIT_DIR / "pipeline.joblib"
        or model_sha256 != MODEL_SHA256
        or (root / AUTH_PATH).resolve(strict=True) != root / AUTH_PATH
        or (root / CLAIM_PATH).resolve(strict=True) != root / CLAIM_PATH
    ):
        raise AuthorizationError("Real test paths or model identity differ")
    verify_approved_fit(root)
    if (
        config != yaml.safe_load((root / CONFIG_PATH).read_bytes())
        or fit_manifest != json.loads((root / FIT_DIR / "training_manifest.json").read_bytes())
        or subprocess.check_output(["git", "status", "--porcelain"], cwd=root, text=True).strip()
    ):
        raise AuthorizationError("Claimed test inputs or checkout changed")
    raw = (root / AUTH_PATH).read_bytes()
    if hashlib.sha256(raw).hexdigest() != authorization_sha256:
        raise AuthorizationError("Claimed test authorization checksum differs")
    record = json.loads(raw)
    expected = expected_test_authorization(
        root,
        source_revision,
        sha256_file(root / CONFIG_PATH),
        sha256_file(root / "data/processed/split_v1.json"),
        record.get("attempt_id"),
        record.get("approval_reference"),
    )
    if (
        record != expected
        or subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root, text=True).strip()
        != source_revision
    ):
        raise AuthorizationError("Claimed test approval differs from current revision")
    _verify_existing_claim(root / CLAIM_PATH, _claim_payload(expected, authorization_sha256))
    return expected


def execute_test_only(
    source: Source,
    config: dict[str, Any],
    output: Path,
    pipeline_path: Path,
    fit_manifest: dict[str, Any],
    *,
    authorization_sha256: str,
    model_sha256: str,
    source_revision: str,
    attempt_id: str | None = None,
    synthetic_model_loader: Callable[[Path], Any] | None = None,
) -> dict[str, Any]:
    """One ordered prediction and score. Caller must have claimed the attempt."""
    if type(source) is RealSource:
        if synthetic_model_loader is not None:
            raise AuthorizationError("Production cannot accept an injected model loader")
        approved = _verify_claimed_real_context(
            source,
            output,
            pipeline_path,
            fit_manifest,
            config,
            authorization_sha256,
            model_sha256,
            source_revision,
        )
        if attempt_id != approved["attempt_id"]:
            raise AuthorizationError("Production attempt ID differs from approved claim")
    elif type(source) not in (SyntheticFixtureSource, SyntheticFileSource):
        raise AuthorizationError("Only the fixed real source or synthetic fixture is accepted")
    if attempt_id is None:
        if type(source) is RealSource:
            raise AuthorizationError("Production test attempt ID is required")
        attempt_id = "synthetic-attempt"
    state = RunState(
        output,
        attempt_id=attempt_id,
        authorization_sha256=authorization_sha256,
        source_revision=source_revision,
        model_sha256=model_sha256,
    )  # Exclusive directory creation is the second one-shot barrier.
    stage = "reserved"
    try:
        snapshot = source.preflight(config)
        if (
            len(snapshot.fit_pairs) != config["partitions"]["final_fit_count"]
            or len(snapshot.test_pairs) != config["partitions"]["final_test_count"]
        ):
            raise ValueError("Frozen population count differs")
        state.transition("test_access_started")
        stage = "test_access_started"
        source.begin_test_access()
        snapshot = source.verify_full_provenance(snapshot)
        if sha256_file(pipeline_path) != model_sha256:
            raise ValueError("Pipeline changed before deserialization")
        pipeline = (
            joblib.load(pipeline_path)
            if type(source) is RealSource
            else (synthetic_model_loader or joblib.load)(pipeline_path)
        )
        if sha256_file(pipeline_path) != model_sha256:
            raise ValueError("Pipeline changed during deserialization")
        classifier = _validate_pipeline(pipeline, snapshot, fit_manifest)
        test = source.load_test_expression(snapshot)
        pairs = tuple(zip(test.patient_ids, test.sample_ids, strict=True))
        if (
            pairs != snapshot.test_pairs
            or len(set(test.patient_ids)) != len(pairs)
            or len(set(test.sample_ids)) != len(pairs)
            or tuple(test.X.index) != test.sample_ids
            or tuple(test.X.columns) != snapshot.gene_ids
            or not np.isfinite(test.X.to_numpy()).all()
        ):
            raise ValueError("Held-out identities, gene order or values differ")
        predicted = pipeline.predict(test.X)
        raw_probabilities = pipeline.predict_proba(test.X)
        positions = [list(classifier.classes_).index(label) for label in LABELS]
        probabilities = raw_probabilities[:, positions]
        count = len(pairs)
        if (
            len(predicted) != count
            or probabilities.shape != (count, len(LABELS))
            or not set(predicted).issubset(set(LABELS))
            or not np.isfinite(probabilities).all()
            or np.any(probabilities < 0)
            or np.any(probabilities > 1)
            or not np.allclose(probabilities.sum(axis=1), 1, atol=1e-9, rtol=0)
        ):
            raise ValueError("Predictions or aligned probabilities are invalid")
        rows = [
            {
                "patient_id": test.patient_ids[i],
                "sample_id": test.sample_ids[i],
                "predicted_label": str(predicted[i]),
                **{
                    f"probability_{label}": float(probabilities[i, j])
                    for j, label in enumerate(LABELS)
                },
            }
            for i in range(count)
        ]
        prediction_path = output / "test_predictions.tsv"
        _write_predictions(prediction_path, rows, LABELS)
        prediction_sha = _verify_saved_predictions(prediction_path, rows, LABELS)
        test_identity_sha = identity_digest(snapshot.test_pairs)
        seal = {
            "status": "sealed_before_labels",
            "attempt_id": attempt_id,
            "authorization_sha256": authorization_sha256,
            "source_revision": source_revision,
            "model_sha256": model_sha256,
            "test_identity_sha256": test_identity_sha,
            "prediction_sha256": prediction_sha,
            "row_count": count,
            "class_order": list(LABELS),
        }
        _atomic_write_json_new(output / "prediction_seal.json", seal)
        if json.loads((output / "prediction_seal.json").read_bytes()) != seal:
            raise ValueError("Prediction seal differs")
        state.transition(
            "predictions_sealed",
            prediction_sha256=prediction_sha,
            seal_sha256=sha256_file(output / "prediction_seal.json"),
            test_identity_sha256=test_identity_sha,
        )
        stage = "predictions_sealed"
        source.note_predictions_saved(output, snapshot)
        if (
            _verify_saved_predictions(prediction_path, rows, LABELS) != prediction_sha
            or sha256_file(output / "prediction_seal.json")
            != json.loads((output / "run_state.jsonl").read_text().splitlines()[-1])["seal_sha256"]
        ):
            raise ValueError("Predictions changed before label access")
        snapshot = source.verify_label_provenance(snapshot)
        if sha256_file(prediction_path) != prediction_sha:
            raise ValueError("Predictions changed during label provenance verification")
        y_true = source.load_test_labels(snapshot)
        if len(y_true) != count or not set(y_true).issubset(set(LABELS)):
            raise ValueError("Held-out label schema differs")
        if _verify_saved_predictions(prediction_path, rows, LABELS) != prediction_sha:
            raise ValueError("Predictions changed after sealing")
        metrics = _score(y_true, predicted, LABELS)
        if type(source) is RealSource and metrics["per_class"]["Normal-like"]["support"] != 5:
            raise ValueError("Normal-like support differs from frozen test split")
        _write_json_new(output / "test_metrics.json", metrics)
        _write_json_new(
            output / "confusion_matrices.json",
            {
                "class_order": list(LABELS),
                "rows": "true_label",
                "columns": "predicted_label",
                "raw": metrics["confusion_matrix"],
                "row_normalized_fraction": metrics["row_normalized_confusion_matrix"],
            },
        )
        names = (
            "test_predictions.tsv",
            "prediction_seal.json",
            "test_metrics.json",
            "confusion_matrices.json",
        )
        manifest = {
            "version": "final_test_only_v1",
            "scope": "one_shot_final_test",
            "source_revision": source_revision,
            "attempt_id": attempt_id,
            "authorization_sha256": authorization_sha256,
            "fit_revision": FIT_REVISION,
            "fit_attempt_id": FIT_ATTEMPT,
            "fit_authorization_sha256": FIT_AUTH_SHA256,
            "fit_manifest_sha256": MANIFEST_SHA256,
            "model_sha256": model_sha256,
            "sealed_prediction_sha256": prediction_sha,
            "test_identity_sha256": test_identity_sha,
            "input_hashes": snapshot.input_hashes,
            "class_order": list(LABELS),
            "test_count": count,
            "normal_like_small_support_note": (
                "Five documented Normal-like test cases; one case changes recall "
                "by 20 percentage points."
            ),
            "output_hashes": {name: sha256_file(output / name) for name in names},
        }
        _atomic_write_json_new(output / "manifest.json", manifest)
        _verify_output_hashes(output, manifest["output_hashes"])
        if sha256_file(prediction_path) != prediction_sha:
            raise ValueError("Sealed predictions changed before completion")
        manifest_sha = sha256_file(output / "manifest.json")
        _atomic_write_json_new(
            output / "success.json",
            {
                "status": "sealed",
                "manifest_sha256": manifest_sha,
            },
        )
        state.transition(
            "completed",
            manifest_sha256=manifest_sha,
            success_sha256=sha256_file(output / "success.json"),
        )
        verify_completed_test_only(output)
        return manifest
    except Exception as exc:
        _write_json_new(
            output / "failure.json",
            {
                "stage": stage,
                "error_type": type(exc).__name__,
                "message": str(exc),
            },
        )
        state.transition("failed", failed_stage=stage, error_type=type(exc).__name__)
        raise


def run_real(root: Path, authorization_path: Path) -> dict[str, Any]:
    """Future authorized entry point; never creates authorization or retries a claim."""
    root = root.resolve(strict=True)
    if (root / OUTPUT_DIR).exists() or (root / CLAIM_PATH).exists():
        raise FileExistsError("Test attempt or output already exists")
    verify_approved_fit(root)
    record, auth_sha = verify_test_authorization(root, authorization_path)
    source = RealSource(root)
    config = yaml.safe_load((root / CONFIG_PATH).read_bytes())
    validate_protocol(config)
    marker = root / CLAIM_PATH
    claim = {
        "attempt_id": record["attempt_id"],
        "approval_reference": record["approval_reference"],
        "approved_commit": record["approved_commit"],
        "authorization_sha256": auth_sha,
        "canonical_checkout": record["canonical_checkout"],
        "canonical_output": record["canonical_output"],
        "status": "claimed_before_test_access",
    }
    _create_local_claim(marker, claim)
    _verify_existing_claim(marker, _claim_payload(record, auth_sha))
    if sha256_file(root / AUTH_PATH) != auth_sha:
        raise AuthorizationError("Test authorization changed after claim")
    fit_manifest = json.loads((root / FIT_DIR / "training_manifest.json").read_bytes())
    return execute_test_only(
        source,
        config,
        root / OUTPUT_DIR,
        root / FIT_DIR / "pipeline.joblib",
        fit_manifest,
        authorization_sha256=auth_sha,
        model_sha256=MODEL_SHA256,
        source_revision=record["approved_commit"],
        attempt_id=record["attempt_id"],
    )
