"""One-shot final-test infrastructure; real execution requires a separate approval record."""

from __future__ import annotations

import csv
import hashlib
import io
import json
import math
import os
import platform
import stat
import subprocess
import time
import warnings
from collections import Counter
from dataclasses import dataclass, replace
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Protocol

import joblib
import numpy as np
import pandas as pd
import sklearn
from sklearn.exceptions import ConvergenceWarning

from oncorna.final_fit_export import (
    ArrayHeader,
    ExportPlan,
    inspect_npy_header_handle,
    verify_completed_export,
)
from oncorna.modeling import make_pipeline, scores
from oncorna.preprocessing import (
    identity_digest,
    sha256_file,
    verify_historical_config,
    verify_split_provenance,
)
from oncorna.synthetic_paths import (
    synthetic_input_path,
    verify_synthetic_handle,
    verify_synthetic_input,
)

LABELS = ("Luminal A", "Luminal B", "Basal-like", "HER2-enriched", "Normal-like")
FIT_COUNTS = {
    "Luminal A": 337,
    "Luminal B": 153,
    "Basal-like": 113,
    "HER2-enriched": 54,
    "Normal-like": 18,
}
STATES = (
    "reserved",
    "fitting",
    "fitted",
    "test_access_started",
    "predictions_sealed",
    "completed",
    "failed",
)
TEST_CONTAINING_INPUTS = frozenset({"all_gene_expression", "source_expression", "cohort"})
FIT_EXPRESSION_PATH = Path("data/processed/final_fit_v1/fit_expression.npy")
FIT_EXPRESSION_MANIFEST_PATH = Path("data/processed/final_fit_v1/manifest.json")


class AuthorizationError(RuntimeError):
    """A real final-test execution was not explicitly authorized."""


@dataclass(frozen=True)
class FitData:
    X: pd.DataFrame
    y: np.ndarray
    patient_ids: tuple[str, ...]
    sample_ids: tuple[str, ...]


@dataclass(frozen=True)
class TestData:
    X: pd.DataFrame
    patient_ids: tuple[str, ...]
    sample_ids: tuple[str, ...]


@dataclass(frozen=True)
class Snapshot:
    fit_pairs: tuple[tuple[str, str], ...]
    test_pairs: tuple[tuple[str, str], ...]
    gene_ids: tuple[str, ...]
    sample_ids: tuple[str, ...]
    input_hashes: dict[str, str]
    historical_config_compatibility: str | None = None


class Source(Protocol):
    is_synthetic: bool

    def preflight(self, config: dict[str, Any]) -> Snapshot: ...
    def load_fit(self, snapshot: Snapshot) -> FitData: ...
    def load_test_expression(self, snapshot: Snapshot) -> TestData: ...
    def load_test_labels(self, snapshot: Snapshot) -> np.ndarray: ...
    def note_predictions_saved(
        self, output: Path | None = None, snapshot: Snapshot | None = None
    ) -> None: ...
    def begin_test_access(self) -> None: ...
    def verify_full_provenance(self, snapshot: Snapshot) -> Snapshot: ...
    def verify_label_provenance(self, snapshot: Snapshot) -> Snapshot: ...


class SyntheticFixtureSource:
    """In-memory synthetic-only source; no filesystem input paths or real IDs."""

    is_synthetic = True

    def __init__(
        self,
        fit: FitData,
        test: TestData,
        test_labels: np.ndarray,
        *,
        fail_stage: str | None = None,
    ) -> None:
        all_ids = fit.patient_ids + fit.sample_ids + test.patient_ids + test.sample_ids
        if not all(value.startswith("SYN-") for value in all_ids):
            raise ValueError("Synthetic fixture IDs must begin with SYN-")
        self.fit = fit
        self.test = test
        self.test_labels = test_labels
        self.fail_stage = fail_stage
        self.events: list[str] = []

    def preflight(self, config: dict[str, Any]) -> Snapshot:
        self.events.append("preflight")
        if self.fail_stage == "preflight":
            raise RuntimeError("Synthetic preflight failure")
        pairs = tuple(zip(self.fit.patient_ids, self.fit.sample_ids, strict=True))
        test_pairs = tuple(zip(self.test.patient_ids, self.test.sample_ids, strict=True))
        if tuple(self.fit.X.columns) != tuple(self.test.X.columns):
            raise ValueError("Synthetic feature schemas differ")
        return Snapshot(
            pairs,
            test_pairs,
            tuple(self.fit.X.columns),
            self.fit.sample_ids + self.test.sample_ids,
            {"synthetic_fixture": hashlib.sha256(self.fit.X.to_numpy().tobytes()).hexdigest()},
        )

    def load_fit(self, snapshot: Snapshot) -> FitData:
        self.events.append("fit_loaded")
        if self.fail_stage == "fit_load":
            raise RuntimeError("Synthetic fit-load failure")
        return self.fit

    def load_test_expression(self, snapshot: Snapshot) -> TestData:
        self.events.append("test_expression_loaded")
        if self.fail_stage == "test_expression":
            raise RuntimeError("Synthetic test-expression failure")
        return self.test

    def note_predictions_saved(
        self, output: Path | None = None, snapshot: Snapshot | None = None
    ) -> None:
        self.events.append("predictions_saved")

    def begin_test_access(self) -> None:
        self.events.append("test_access_authorized")

    def verify_full_provenance(self, snapshot: Snapshot) -> Snapshot:
        self.events.append("full_provenance_verified")
        return snapshot

    def verify_label_provenance(self, snapshot: Snapshot) -> Snapshot:
        if "predictions_saved" not in self.events:
            raise RuntimeError("Synthetic label provenance requires sealed predictions")
        self.events.append("label_provenance_verified")
        return snapshot

    def load_test_labels(self, snapshot: Snapshot) -> np.ndarray:
        self.events.append("test_labels_loaded")
        if "predictions_saved" not in self.events or self.events[-2] not in {
            "predictions_saved",
            "label_provenance_verified",
        }:
            raise RuntimeError("Synthetic labels accessed before predictions were saved")
        if self.fail_stage == "test_labels":
            raise RuntimeError("Synthetic test-label failure")
        return self.test_labels


class SyntheticFileSource(SyntheticFixtureSource):
    """Synthetic label-file fixture with observable, seal-gated file opens."""

    def __init__(self, fit: FitData, test: TestData, label_path: Path) -> None:
        super().__init__(fit, test, np.asarray([], dtype=object))
        self.label_path = label_path
        self.label_identity = synthetic_input_path(label_path)

    def note_predictions_saved(
        self, output: Path | None = None, snapshot: Snapshot | None = None
    ) -> None:
        if output is None or snapshot is None:
            raise RuntimeError("Synthetic file labels require a saved prediction seal")
        entries = [
            json.loads(line) for line in (output / "run_state.jsonl").read_text().splitlines()
        ]
        seal = json.loads((output / "prediction_seal.json").read_bytes())
        if (
            entries[-1]["state"] != "predictions_sealed"
            or entries[-1]["seal_sha256"] != sha256_file(output / "prediction_seal.json")
            or seal["prediction_sha256"] != sha256_file(output / "test_predictions.tsv")
            or seal["test_identity_sha256"] != identity_digest(snapshot.test_pairs)
            or seal["attempt_id"] != entries[-1]["attempt_id"]
        ):
            raise RuntimeError("Synthetic prediction seal differs")
        self.events.append("predictions_sealed_verified")
        super().note_predictions_saved(output, snapshot)

    def load_test_labels(self, snapshot: Snapshot) -> np.ndarray:
        if (
            "predictions_sealed_verified" not in self.events
            or "label_provenance_verified" not in self.events
        ):
            raise RuntimeError("Synthetic file labels accessed before predictions were sealed")
        verify_synthetic_input(self.label_path, self.label_identity)
        with self.label_path.open("r", encoding="utf-8", newline="") as handle:
            verify_synthetic_handle(handle, self.label_identity)
            verify_synthetic_input(self.label_path, self.label_identity)
            self.events.append("label_file_opened")
            rows = list(csv.DictReader(handle, delimiter="\t"))
        pairs = tuple((row["patient_id"], row["sample_id"]) for row in rows)
        if pairs != snapshot.test_pairs:
            raise ValueError("Synthetic label file identities differ")
        self.events.append("test_labels_loaded")
        return np.asarray([row["true_label"] for row in rows], dtype=object)


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sequence_digest(values: tuple[str, ...]) -> str:
    payload = json.dumps(values, ensure_ascii=False, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def validate_memberships(
    partitions: dict[str, tuple[tuple[str, str], ...]],
    expected: dict[str, int],
) -> tuple[tuple[tuple[str, str], ...], tuple[tuple[str, str], ...]]:
    if set(partitions) != {"train", "validation", "test"}:
        raise ValueError("Frozen split must have exactly train, validation, and test")
    seen_patients: set[str] = set()
    seen_samples: set[str] = set()
    for name in ("train", "validation", "test"):
        pairs = partitions[name]
        if len(pairs) != expected[name]:
            raise ValueError(f"Frozen {name} count differs from approved protocol")
        for patient, sample in pairs:
            if not patient or not sample or patient in seen_patients or sample in seen_samples:
                raise ValueError("Frozen split has missing or overlapping patient/sample IDs")
            seen_patients.add(patient)
            seen_samples.add(sample)
    fit_pairs = partitions["train"] + partitions["validation"]
    if len(fit_pairs) != expected["train"] + expected["validation"]:
        raise ValueError("Final fitting membership is incomplete")
    return fit_pairs, partitions["test"]


def validate_protocol(config: dict[str, Any], *, synthetic_fixture: bool = False) -> None:
    p = config["pipeline"]
    r = config["reporting"]
    if (
        config["version"] != "final_test_v1"
        or config["execution_requires_separate_human_approval"] is not True
        or config["partitions"]["final_fit"] != ["train", "validation"]
        or config["partitions"]["final_fit_count"]
        != (config["partitions"]["original_train"] + config["partitions"]["original_validation"])
        or config["partitions"]["final_test_count"] < 1
        or config["partitions"]["split_membership_must_remain_unchanged"] is not True
        or tuple(config["task"]["class_order"]) != LABELS
        or tuple(r["class_order"]) != LABELS
        or p["gene_filter"]["expression_gt"] != 1.0
        or p["gene_filter"]["minimum_sample_fraction"] != 0.20
        or p["gene_filter"]["required_fit_sample_count"]
        != math.ceil(0.20 * config["partitions"]["final_fit_count"])
        or p["gene_filter"]["refit_on_final_fit_patients"] is not True
        or p["scaler"] != "StandardScaler"
        or any(p[key] is not None for key in ("imputer", "feature_selector", "pca", "calibration"))
        or p["additional_hyperparameter_search"] is not False
        or p["classifier"]
        != {
            "type": "multinomial_logistic_regression",
            "penalty": "l2",
            "C": 1.0,
            "solver": "lbfgs",
            "fit_intercept": True,
            "class_weight": None,
            "tol": 0.0001,
            "max_iter": 1000,
        }
        or config["metrics"]["primary"] != "macro_f1"
        or config["metrics"]["zero_division"] != 0
        or r["zero_division"] != 0
        or r["confusion_rows"] != "true_label"
        or r["confusion_columns"] != "predicted_label"
        or r["row_normalization"] != "divide_each_row_by_its_true_class_support"
        or r["row_normalized_machine_units"] != "fraction_0_to_1"
        or r["row_normalized_human_units"] != "percent"
        or r["machine_numeric_precision"] != "full_float64_json"
        or r["human_metric_decimal_places"] != 4
        or r["human_percentage_decimal_places"] != 1
        or r["probability_columns"] != "ordered_by_reporting_class_order"
        or r["probability_row_sum_absolute_tolerance"] != 1.0e-9
        or r["support_reconciliation"] != "rows_equal_true_class_support_and_total_test_count"
        or r["prediction_reconciliation"]
        != "columns_equal_predicted_class_counts_and_total_test_count"
        or config["outputs"]["exclusive_create"] is not True
        or config["outputs"]["overwrite"] is not False
        or config["outputs"]["patient_level_and_model_artifacts"] != "local_only"
    ):
        raise ValueError("Final-test configuration differs from approved scientific protocol")
    if not synthetic_fixture and (
        config["partitions"]["original_train"] != 506
        or config["partitions"]["original_validation"] != 169
        or config["partitions"]["final_fit_count"] != 675
        or config["partitions"]["final_test_count"] != 169
        or config["task"]["feature_count_before_filter"] != 20530
    ):
        raise ValueError("Real final-test population differs from frozen protocol")


def modeling_parameters(config: dict[str, Any]) -> dict[str, Any]:
    classifier = config["pipeline"]["classifier"]
    return {
        "expression_gt": config["pipeline"]["gene_filter"]["expression_gt"],
        "minimum_sample_fraction": config["pipeline"]["gene_filter"]["minimum_sample_fraction"],
        "penalty": classifier["penalty"],
        "solver": classifier["solver"],
        "class_weight": classifier["class_weight"],
        "fit_intercept": classifier["fit_intercept"],
        "max_iter": classifier["max_iter"],
        "tol": classifier["tol"],
    }


def _write_json_new(path: Path, value: Any) -> None:
    with path.open("x", encoding="utf-8", newline="\n") as handle:
        json.dump(value, handle, ensure_ascii=False, indent=2, sort_keys=True)
        handle.write("\n")


def _atomic_write_json_new(path: Path, value: Any) -> None:
    if path.exists():
        raise FileExistsError(f"Sealed output already exists: {path}")
    temporary = path.with_name(f".{path.name}.pending")
    with temporary.open("x", encoding="utf-8", newline="\n") as handle:
        json.dump(value, handle, ensure_ascii=False, indent=2, sort_keys=True)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def _verify_output_hashes(output: Path, expected: dict[str, str]) -> None:
    for name, digest in expected.items():
        if sha256_file(output / name) != digest:
            raise ValueError(f"Sealed output checksum mismatch: {name}")


def verify_completed_run(output: Path) -> dict[str, Any]:
    """Accept success only when the terminal audit state and all seals agree."""
    entries = [
        json.loads(line)
        for line in (output / "run_state.jsonl").read_text(encoding="utf-8").splitlines()
    ]
    if not entries or entries[-1]["state"] != "completed":
        raise ValueError("Final-test attempt has no terminal completed state")
    terminal = entries[-1]
    manifest_path = output / "manifest.json"
    success_path = output / "success.json"
    if sha256_file(manifest_path) != terminal["manifest_sha256"]:
        raise ValueError("Completed manifest checksum does not match audit state")
    if sha256_file(success_path) != terminal["success_sha256"]:
        raise ValueError("Success marker checksum does not match audit state")
    success = json.loads(success_path.read_text(encoding="utf-8"))
    if success != {"manifest_sha256": terminal["manifest_sha256"], "status": "sealed"}:
        raise ValueError("Success marker does not identify the completed manifest")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    _verify_output_hashes(output, manifest["output_hashes"])
    return manifest


class RunState:
    def __init__(self, output: Path, **binding: Any) -> None:
        output.mkdir(parents=True, exist_ok=False)
        self.path = output / "run_state.jsonl"
        self.sequence = 0
        self.binding = binding
        self.transition("reserved")

    def transition(self, state: str, **details: Any) -> None:
        if state not in STATES:
            raise ValueError(f"Unknown run state: {state}")
        self.sequence += 1
        entry = {
            "sequence": self.sequence,
            "state": state,
            "at_utc": utc_now(),
            **self.binding,
            **details,
        }
        mode = "x" if self.sequence == 1 else "a"
        with self.path.open(mode, encoding="utf-8", newline="\n") as handle:
            handle.write(json.dumps(entry, ensure_ascii=False, sort_keys=True) + "\n")
            handle.flush()
            os.fsync(handle.fileno())


def _write_predictions(path: Path, rows: list[dict[str, Any]], labels: tuple[str, ...]) -> None:
    fields = ["patient_id", "sample_id", "predicted_label", *(f"probability_{x}" for x in labels)]
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
        handle.flush()
        os.fsync(handle.fileno())


def _verify_saved_predictions(
    path: Path, rows: list[dict[str, Any]], labels: tuple[str, ...]
) -> str:
    with path.open(encoding="utf-8", newline="") as handle:
        saved = list(csv.DictReader(handle, delimiter="\t"))
    if len(saved) != len(rows):
        raise ValueError("Saved prediction count differs from frozen test count")
    fields = ["patient_id", "sample_id", "predicted_label", *(f"probability_{x}" for x in labels)]
    if not saved or list(saved[0]) != fields:
        raise ValueError("Saved prediction column order differs from frozen class order")
    if len({row["patient_id"] for row in saved}) != len(rows) or len(
        {row["sample_id"] for row in saved}
    ) != len(rows):
        raise ValueError("Saved predictions contain duplicate patient or sample IDs")
    for actual, expected in zip(saved, rows, strict=True):
        if any(actual[key] != str(expected[key]) for key in fields):
            raise ValueError("Saved predictions differ from in-memory predictions")
    return sha256_file(path)


def _score(y_true: np.ndarray, y_pred: np.ndarray, labels: tuple[str, ...]) -> dict[str, Any]:
    result = scores(y_true, y_pred, list(labels))
    matrix = np.asarray(result["confusion_matrix"], dtype=np.int64)
    support = [result["per_class"][label]["support"] for label in labels]
    predicted_counts = {label: int(np.count_nonzero(y_pred == label)) for label in labels}
    if (
        matrix.shape != (5, 5)
        or matrix.sum(axis=1).tolist() != support
        or matrix.sum(axis=0).tolist() != list(predicted_counts.values())
        or int(matrix.sum()) != len(y_true)
    ):
        raise ValueError("Class support or predicted counts do not reconcile")
    if any(value == 0 for value in support):
        raise ValueError("Frozen test lacks a required class")
    result["row_normalized_confusion_matrix"] = (
        matrix / matrix.sum(axis=1, keepdims=True)
    ).tolist()
    result["predicted_class_counts"] = predicted_counts
    result["class_order"] = list(labels)
    return result


def execute_once(
    source: Source,
    config: dict[str, Any],
    output: Path,
    *,
    config_sha256: str,
    source_revision: str,
    authorization_sha256: str | None,
) -> dict[str, Any]:
    """Execute stages once; only a separately authorized caller may pass a real source."""
    synthetic = type(source) is SyntheticFixtureSource
    validate_protocol(config, synthetic_fixture=synthetic)
    if not synthetic:
        raise AuthorizationError("Real final-test execution is disabled in Phase 2D-B")
    try:
        snapshot = source.preflight(config)
        expected_fit = config["partitions"]["final_fit_count"]
        expected_test = config["partitions"]["final_test_count"]
        fit_pairs = snapshot.fit_pairs
        test_pairs = snapshot.test_pairs
        if (
            len(fit_pairs) != expected_fit
            or len(test_pairs) != expected_test
            or set(fit_pairs) & set(test_pairs)
            or len({p for p, _ in fit_pairs + test_pairs}) != expected_fit + expected_test
            or len({s for _, s in fit_pairs + test_pairs}) != expected_fit + expected_test
        ):
            raise ValueError("Preflight partition membership failed")
    except Exception as exc:
        state = RunState(output)
        _write_json_new(
            output / "failure.json",
            {
                "stage": "preflight",
                "error_type": type(exc).__name__,
                "message": str(exc),
                "at_utc": utc_now(),
                "source_revision": source_revision,
                "config_sha256": config_sha256,
            },
        )
        state.transition("failed", failed_stage="preflight", error_type=type(exc).__name__)
        raise
    state = RunState(output)
    started = time.perf_counter()
    stage = "reserved"
    try:
        state.transition("fitting")
        stage = "fitting"
        fit = source.load_fit(snapshot)
        if (
            tuple(zip(fit.patient_ids, fit.sample_ids, strict=True)) != fit_pairs
            or tuple(fit.X.index) != fit.sample_ids
            or tuple(fit.X.columns) != snapshot.gene_ids
            or len(fit.y) != expected_fit
            or not np.isfinite(fit.X.to_numpy()).all()
        ):
            raise ValueError("Final fitting data differ from frozen membership or gene schema")
        if not synthetic and Counter(fit.y) != Counter(FIT_COUNTS):
            raise ValueError("Final fitting class counts differ from frozen split")
        pipeline = make_pipeline(modeling_parameters(config), 1.0)
        fit_started = time.perf_counter()
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            pipeline.fit(fit.X, fit.y)
        fit_seconds = time.perf_counter() - fit_started
        gene_filter = pipeline.named_steps["gene_filter"]
        scaler = pipeline.named_steps["scaler"]
        classifier = pipeline.named_steps["classifier"]
        warning_rows = [
            {"category": w.category.__name__, "message": str(w.message)} for w in caught
        ]
        if any(issubclass(w.category, ConvergenceWarning) for w in caught):
            raise RuntimeError("ConvergenceWarning: final fit stopped before test access")
        required = math.ceil(
            config["pipeline"]["gene_filter"]["minimum_sample_fraction"] * expected_fit
        )
        if (
            gene_filter.fitting_sample_ids_ != fit.sample_ids
            or gene_filter.required_count_ != required
            or np.any(np.asarray(scaler.n_samples_seen_) != expected_fit)
            or tuple(classifier.classes_) != tuple(sorted(LABELS))
            or not gene_filter.selected_gene_ids_
        ):
            raise ValueError("Fitted pipeline failed training-only verification")
        gene_ids = tuple(gene_filter.selected_gene_ids_)
        joblib.dump(pipeline, output / "pipeline.joblib")
        with (output / "training_gene_list.tsv").open("x", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["gene_id"])
            writer.writerows((gene,) for gene in gene_ids)
        fit_manifest = {
            "fit_count": expected_fit,
            "fit_identity_sha256": identity_digest(fit_pairs),
            "retained_gene_count": len(gene_ids),
            "ordered_gene_list_sha256": sequence_digest(gene_ids),
            "required_prevalence_count": required,
            "scaler_fit_count": expected_fit,
            "class_order": list(LABELS),
            "solver_iterations": classifier.n_iter_.tolist(),
            "warnings": warning_rows,
            "fit_seconds": fit_seconds,
        }
        _write_json_new(output / "training_manifest.json", fit_manifest)
        state.transition("fitted", gene_count=len(gene_ids))
        stage = "fitted"
        # This stage remains a separate gate even after the earlier real-run
        # block is removed in a later reviewed revision.
        if not synthetic:
            raise AuthorizationError("Real test-access gate remains disabled")
        state.transition("test_access_started")
        stage = "test_access_started"
        source.begin_test_access()
        snapshot = source.verify_full_provenance(snapshot)
        test = source.load_test_expression(snapshot)
        if (
            tuple(zip(test.patient_ids, test.sample_ids, strict=True)) != test_pairs
            or tuple(test.X.index) != test.sample_ids
            or tuple(test.X.columns) != snapshot.gene_ids
            or not np.isfinite(test.X.to_numpy()).all()
        ):
            raise ValueError("Test expression differs from frozen membership or gene schema")
        predicted = pipeline.predict(test.X)
        raw_probabilities = pipeline.predict_proba(test.X)
        positions = [list(classifier.classes_).index(label) for label in LABELS]
        probabilities = raw_probabilities[:, positions]
        tolerance = config["reporting"]["probability_row_sum_absolute_tolerance"]
        if (
            len(predicted) != expected_test
            or probabilities.shape != (expected_test, len(LABELS))
            or not set(predicted).issubset(set(LABELS))
            or not np.isfinite(probabilities).all()
            or np.any(probabilities < 0)
            or np.any(probabilities > 1)
            or not np.allclose(probabilities.sum(axis=1), 1, atol=tolerance, rtol=0)
        ):
            raise ValueError("Test predictions or aligned probabilities are invalid")
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
            for i in range(expected_test)
        ]
        _write_predictions(output / "test_predictions.tsv", rows, LABELS)
        prediction_sha256 = _verify_saved_predictions(output / "test_predictions.tsv", rows, LABELS)
        source.note_predictions_saved()
        y_true = source.load_test_labels(snapshot)
        if len(y_true) != expected_test or not set(y_true).issubset(set(LABELS)):
            raise ValueError("Test labels differ from frozen class schema")
        metric = _score(y_true, predicted, LABELS)
        _write_json_new(output / "test_metrics.json", metric)
        _write_json_new(
            output / "confusion_matrices.json",
            {
                "class_order": list(LABELS),
                "rows": "true_label",
                "columns": "predicted_label",
                "raw": metric["confusion_matrix"],
                "row_normalized_fraction": metric["row_normalized_confusion_matrix"],
            },
        )
        output_names = (
            "pipeline.joblib",
            "training_gene_list.tsv",
            "training_manifest.json",
            "test_predictions.tsv",
            "test_metrics.json",
            "confusion_matrices.json",
        )
        manifest = {
            "version": "final_test_v1",
            "scope": "one_shot_final_test",
            "source_revision": source_revision,
            "config_sha256": config_sha256,
            "authorization_sha256": authorization_sha256,
            "input_hashes": snapshot.input_hashes,
            "historical_config_compatibility": snapshot.historical_config_compatibility,
            "frozen_split_sha256": snapshot.input_hashes.get("frozen_split"),
            "fit_identity_sha256": fit_manifest["fit_identity_sha256"],
            "ordered_gene_list_sha256": fit_manifest["ordered_gene_list_sha256"],
            "model_sha256": sha256_file(output / "pipeline.joblib"),
            "sealed_prediction_sha256": prediction_sha256,
            "software": {
                "python": platform.python_version(),
                "numpy": np.__version__,
                "pandas": pd.__version__,
                "scikit_learn": sklearn.__version__,
                "joblib": joblib.__version__,
            },
            "solver_iterations": classifier.n_iter_.tolist(),
            "warnings": warning_rows,
            "fit_seconds": fit_seconds,
            "elapsed_seconds": time.perf_counter() - started,
            "finished_at_utc": utc_now(),
            "output_hashes": {name: sha256_file(output / name) for name in output_names},
        }
        _atomic_write_json_new(output / "manifest.json", manifest)
        if json.loads((output / "manifest.json").read_text(encoding="utf-8")) != manifest:
            raise ValueError("Final manifest differs from sealed in-memory record")
        _verify_output_hashes(output, manifest["output_hashes"])
        manifest_sha256 = sha256_file(output / "manifest.json")
        success = {"manifest_sha256": manifest_sha256, "status": "sealed"}
        _atomic_write_json_new(output / "success.json", success)
        if json.loads((output / "success.json").read_text(encoding="utf-8")) != success:
            raise ValueError("Success marker differs from sealed manifest identity")
        success_sha256 = sha256_file(output / "success.json")
        state.transition(
            "completed", manifest_sha256=manifest_sha256, success_sha256=success_sha256
        )
        verify_completed_run(output)
        return manifest
    except Exception as exc:
        _write_json_new(
            output / "failure.json",
            {
                "stage": stage,
                "error_type": type(exc).__name__,
                "message": str(exc),
                "at_utc": utc_now(),
                "source_revision": source_revision,
                "config_sha256": config_sha256,
            },
        )
        state.transition("failed", failed_stage=stage, error_type=type(exc).__name__)
        raise


def _source_stamp(info: os.stat_result) -> tuple[int, int, int, int, int]:
    return (info.st_dev, info.st_ino, info.st_size, info.st_mtime_ns, info.st_ctime_ns)


def _authenticated_source_rows(
    source_path: Path,
    fit_path: Path,
    snapshot: Snapshot,
    expected_sha256: str,
    recorded_offset: int,
) -> tuple[np.ndarray, tuple[int, int, int, int, int]]:
    """Hash and select fit/test rows from one immutable in-memory copy of source bytes."""
    if source_path.resolve(strict=True) != source_path:
        raise ValueError("Original expression source path is redirected")
    with source_path.open("rb", buffering=0) as handle:
        before = _source_stamp(os.fstat(handle.fileno()))
        if (
            not stat.S_ISREG(os.fstat(handle.fileno()).st_mode)
            or _source_stamp(source_path.stat()) != before
        ):
            raise ValueError("Original source identity differs before authentication")
        header = inspect_npy_header_handle(
            handle, (len(snapshot.sample_ids), len(snapshot.gene_ids))
        )
        if header.data_offset != recorded_offset:
            raise ValueError("Original source array offset differs from export")
        handle.seek(0)
        raw = handle.read()
        if (
            _source_stamp(os.fstat(handle.fileno())) != before
            or _source_stamp(source_path.stat()) != before
            or hashlib.sha256(raw).hexdigest() != expected_sha256
        ):
            raise ValueError("Original source changed or checksum differs")
    values = np.load(io.BytesIO(raw), allow_pickle=False)
    if values.shape != header.shape or values.dtype != np.dtype("<f8"):
        raise ValueError("Authenticated source schema differs")
    fit_values = np.load(fit_path, mmap_mode="r", allow_pickle=False)
    if fit_values.shape != (len(snapshot.fit_pairs), len(snapshot.gene_ids)):
        raise ValueError("Fitting-only source schema differs")
    positions = {sample: i for i, sample in enumerate(snapshot.sample_ids)}
    fit_rows = [positions[sample] for _, sample in snapshot.fit_pairs]
    test_rows = [positions[sample] for _, sample in snapshot.test_pairs]
    if not np.array_equal(np.asarray(fit_values), values[fit_rows, :]):
        raise ValueError("Fitting-only expression differs from authenticated original")
    test_values = np.array(values[test_rows, :], copy=True)
    test_values.setflags(write=False)
    return test_values, before


class RealSource:
    """Lazy file-backed source. Construction never loads expression or labels."""

    is_synthetic = False

    def __init__(self, root: Path) -> None:
        self.root = root
        self.test_access_started = False
        self.full_provenance_verified = False
        self.predictions_saved = False
        self.prediction_sha256: str | None = None
        self.prediction_seal_sha256: str | None = None
        self.authenticated_test_values: np.ndarray | None = None
        self.authenticated_source_stamp: tuple[int, int, int, int, int] | None = None

    def _require_access(self, state: str) -> dict[str, str]:
        from oncorna.final_test_consumer import verify_active_test_attempt

        return verify_active_test_attempt(self.root, state)

    def _path(self, config: dict[str, Any], key: str) -> Path:
        return self.root / config["inputs"][key]

    def preflight(self, config: dict[str, Any]) -> Snapshot:
        self.config = config
        expected = config["inputs"]["expected_sha256"]
        hashes = {}
        for key, digest in expected.items():
            if key in TEST_CONTAINING_INPUTS:
                continue
            path = self._path(config, key)
            actual = sha256_file(path)
            if actual != digest:
                raise ValueError(f"Frozen input checksum mismatch: {key}")
            hashes[key] = actual
        split = json.loads(self._path(config, "frozen_split").read_text(encoding="utf-8"))
        if split.get("split_id") != "split_v1" or split.get("split_unit") != "patient_id":
            raise ValueError("Unexpected frozen split schema")
        parts = {
            name: tuple(
                (x["patient_id"], x["sample_id"])
                for x in split["partitions"][name]["patient_sample_map"]
            )
            for name in ("train", "validation", "test")
        }
        fit_pairs, test_pairs = validate_memberships(
            parts,
            {
                "train": config["partitions"]["original_train"],
                "validation": config["partitions"]["original_validation"],
                "test": config["partitions"]["final_test_count"],
            },
        )
        schema = json.loads(self._path(config, "all_gene_schema").read_text(encoding="utf-8"))
        gene_ids = tuple(schema["gene_ids"])
        sample_ids = tuple(schema["sample_ids"])
        if (
            schema["orientation"] != "samples_by_genes"
            or schema["expression_units"] != config["task"]["expression_units"]
            or schema["shape"]
            != [
                config["partitions"]["final_fit_count"] + config["partitions"]["final_test_count"],
                config["task"]["feature_count_before_filter"],
            ]
            or len(gene_ids) != len(set(gene_ids))
            or len(sample_ids) != len(set(sample_ids))
            or set(sample_ids) != {sample for _, sample in fit_pairs + test_pairs}
        ):
            raise ValueError("All-gene schema differs from frozen population")
        historical = verify_historical_config(
            self.root / "configs/default.yaml", split["configuration"]["sha256"]
        )
        historical_manifest = self.root / "data/processed/ml_matrix_manifest.json"
        expected_historical_manifest = split["expression_matrix_header_validation"][
            "matrix_manifest_sha256"
        ]
        if sha256_file(historical_manifest) != expected_historical_manifest:
            raise ValueError("Historical matrix manifest checksum differs from frozen split")
        hashes["historical_matrix_manifest"] = expected_historical_manifest
        training_predictions = self.root / "data/processed/modeling_v1/oof_predictions.tsv"
        validation_predictions = (
            self.root / "data/processed/validation_v1/validation_predictions.tsv"
        )
        for manifest_key, artifact, hash_key in (
            ("phase2b_manifest", training_predictions, "training_oof_predictions"),
            ("phase2c_manifest", validation_predictions, "validation_predictions"),
        ):
            manifest = json.loads(self._path(config, manifest_key).read_text(encoding="utf-8"))
            recorded = manifest["output_hashes"][artifact.name]
            actual = sha256_file(artifact)
            if actual != recorded:
                raise ValueError(f"Training-side label artifact checksum mismatch: {artifact.name}")
            hashes[hash_key] = actual
        fit_manifest_path = self.root / FIT_EXPRESSION_MANIFEST_PATH
        fit_manifest = json.loads(fit_manifest_path.read_text(encoding="utf-8"))
        fit_expression_path = self.root / FIT_EXPRESSION_PATH
        plan = ExportPlan(
            source_array=self._path(config, "all_gene_expression"),
            output_dir=self.root / FIT_EXPRESSION_PATH.parent,
            fit_pairs=fit_pairs,
            test_pairs=test_pairs,
            sample_ids=sample_ids,
            gene_ids=gene_ids,
            frozen_split_sha256=expected["frozen_split"],
            source_schema_sha256=expected["all_gene_schema"],
            source_all_gene_expression_sha256=expected["all_gene_expression"],
            source_expression_sha256=expected["source_expression"],
            approved_cohort_identifier="TCGA-BRCA Xena primary-tumor cohort",
            approved_cohort_path=config["inputs"]["cohort"],
            approved_cohort_sha256=expected["cohort"],
            expression_units=schema["expression_units"],
            command=fit_manifest.get("command", ""),
        )
        header = ArrayHeader(
            tuple(schema["shape"]),
            fit_manifest.get("source_array_data_offset", -1),
            len(gene_ids) * 8,
        )
        verify_completed_export(fit_expression_path.parent, plan, header)
        fit_array = np.load(fit_expression_path, mmap_mode="r", allow_pickle=False)
        if fit_array.shape != (len(fit_pairs), len(gene_ids)) or fit_array.dtype != np.dtype(
            schema["dtype"]
        ):
            raise ValueError("Fitting-only expression shape or dtype differs from schema")
        hashes["fit_expression"] = fit_manifest["fit_expression_sha256"]
        hashes["fit_expression_manifest"] = sha256_file(fit_manifest_path)
        return Snapshot(
            fit_pairs,
            test_pairs,
            gene_ids,
            sample_ids,
            hashes,
            historical["compatibility_mode"],
        )

    def begin_test_access(self) -> None:
        self._require_access("test_access_started")
        self.test_access_started = True

    def verify_full_provenance(self, snapshot: Snapshot) -> Snapshot:
        """Access full-cohort bytes only after the explicit test-access transition."""
        self._require_access("test_access_started")
        if not self.test_access_started:
            raise AuthorizationError("Full-cohort provenance requires the test-access stage")
        expected = self.config["inputs"]["expected_sha256"]
        hashes = dict(snapshot.input_hashes)
        original_path = self._path(self.config, "source_expression")
        if original_path.resolve(strict=True) != original_path:
            raise AuthorizationError("Original source path is redirected")
        original_sha = sha256_file(original_path)
        if original_sha != expected["source_expression"]:
            raise ValueError("Frozen original-source checksum mismatch")
        hashes["source_expression"] = original_sha
        fit_manifest_path = self.root / FIT_EXPRESSION_MANIFEST_PATH
        if sha256_file(fit_manifest_path) != snapshot.input_hashes["fit_expression_manifest"]:
            raise ValueError("Fitting-only manifest changed after preflight")
        fit_manifest = json.loads(fit_manifest_path.read_text(encoding="utf-8"))
        if sha256_file(self.root / FIT_EXPRESSION_PATH) != snapshot.input_hashes["fit_expression"]:
            raise ValueError("Fitting-only expression changed before source comparison")
        test_values, source_stamp = _authenticated_source_rows(
            self._path(self.config, "all_gene_expression"),
            self.root / FIT_EXPRESSION_PATH,
            snapshot,
            expected["all_gene_expression"],
            fit_manifest["source_array_data_offset"],
        )
        self.authenticated_test_values = test_values
        self.authenticated_source_stamp = source_stamp
        hashes["all_gene_expression"] = expected["all_gene_expression"]
        self.full_provenance_verified = True
        return replace(snapshot, input_hashes=hashes)

    def verify_label_provenance(self, snapshot: Snapshot) -> Snapshot:
        """Open label-containing cohort bytes only after predictions are sealed."""
        self._require_access("predictions_sealed")
        if not self.full_provenance_verified or not self.predictions_saved:
            raise AuthorizationError("Cohort provenance requires sealed predictions")
        self._verify_prediction_seal(snapshot)
        expected = self.config["inputs"]["expected_sha256"]
        cohort_sha = sha256_file(self._path(self.config, "cohort"))
        if cohort_sha != expected["cohort"]:
            raise ValueError("Frozen cohort checksum mismatch")
        split = json.loads(self._path(self.config, "frozen_split").read_text(encoding="utf-8"))
        provenance = verify_split_provenance(
            split,
            cohort_path=self._path(self.config, "cohort"),
            historical_config_path=self.root / "configs/default.yaml",
            historical_matrix_path=self.root / "data/processed/ml_expression_matrix.tsv.gz",
            historical_matrix_manifest_path=self.root / "data/processed/ml_matrix_manifest.json",
        )
        hashes = dict(snapshot.input_hashes)
        hashes["cohort"] = cohort_sha
        hashes["historical_matrix"] = provenance["historical_matrix_sha256"]
        return replace(snapshot, input_hashes=hashes)

    def _fit_expression(self, snapshot: Snapshot) -> pd.DataFrame:
        """Open only the fitting-only artifact; never the full-cohort array."""
        path = self.root / FIT_EXPRESSION_PATH
        if sha256_file(path) != snapshot.input_hashes["fit_expression"]:
            raise ValueError("Fitting-only expression changed after preflight")
        values = np.load(path, mmap_mode="r", allow_pickle=False)
        if values.shape != (len(snapshot.fit_pairs), len(snapshot.gene_ids)):
            raise ValueError("Fitting-only array shape differs from frozen schema")
        selected = np.asarray(values)
        if not np.isfinite(selected).all():
            raise ValueError("Fitting-only expression includes nonfinite values")
        return pd.DataFrame(
            selected,
            index=tuple(sample for _, sample in snapshot.fit_pairs),
            columns=snapshot.gene_ids,
        )

    def _test_expression(self, snapshot: Snapshot) -> pd.DataFrame:
        self._require_access("test_access_started")
        if not self.full_provenance_verified:
            raise AuthorizationError("Test expression requires verified full-cohort provenance")
        path = self._path(self.config, "all_gene_expression")
        if (
            self.authenticated_test_values is None
            or self.authenticated_source_stamp is None
            or path.resolve(strict=True) != path
            or _source_stamp(path.stat()) != self.authenticated_source_stamp
            or sha256_file(path) != snapshot.input_hashes["all_gene_expression"]
        ):
            raise ValueError("Authenticated source changed before test-row use")
        samples = tuple(sample for _, sample in snapshot.test_pairs)
        selected = self.authenticated_test_values
        if selected.shape != (len(samples), len(snapshot.gene_ids)):
            raise ValueError("Authenticated test-row schema differs")
        if not np.isfinite(selected).all():
            raise ValueError("Selected expression includes nonfinite values")
        return pd.DataFrame(selected.copy(), index=samples, columns=snapshot.gene_ids)

    def _fit_labels(self, snapshot: Snapshot) -> np.ndarray:
        """Read only previously sealed training and validation prediction artifacts."""
        sections = (
            (
                self.root / "data/processed/modeling_v1/oof_predictions.tsv",
                snapshot.fit_pairs[: self.config["partitions"]["original_train"]],
                True,
            ),
            (
                self.root / "data/processed/validation_v1/validation_predictions.tsv",
                snapshot.fit_pairs[self.config["partitions"]["original_train"] :],
                False,
            ),
        )
        labels: list[str] = []
        for path, pairs, is_training in sections:
            hash_key = "training_oof_predictions" if is_training else "validation_predictions"
            if sha256_file(path) != snapshot.input_hashes[hash_key]:
                raise ValueError("Fitting-only labels changed after preflight")
            found: dict[tuple[str, str], str] = {}
            with path.open("r", encoding="utf-8", newline="") as handle:
                for row in csv.DictReader(handle, delimiter="\t"):
                    if is_training and (row["model"], row["C"]) != (
                        "logistic_l2_multinomial",
                        "1.0",
                    ):
                        continue
                    pair = (row["patient_id"], row["sample_id"])
                    if pair in found:
                        raise ValueError("Duplicate fitting label in training-side artifact")
                    found[pair] = row["true_label"]
            if set(found) != set(pairs):
                raise ValueError("Training-side label artifact differs from frozen fit membership")
            labels.extend(found[pair] for pair in pairs)
        return np.asarray(labels, dtype=object)

    def _test_labels_from_cohort(self, pairs: tuple[tuple[str, str], ...]) -> np.ndarray:
        self._require_access("predictions_sealed")
        if not self.full_provenance_verified or not self.predictions_saved:
            raise AuthorizationError("Test labels require authorized access and sealed predictions")
        self._verify_prediction_seal(Snapshot((), pairs, (), (), {}))
        wanted = set(pairs)
        found: dict[tuple[str, str], str] = {}
        with self._path(self.config, "cohort").open(
            "r", encoding="utf-8-sig", newline=""
        ) as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                pair = (row["patient_id"], row["sample_id"])
                if pair in wanted:
                    if pair in found:
                        raise ValueError("Duplicate selected cohort patient/sample pair")
                    found[pair] = row["pam50_normalized_label"]
        if set(found) != wanted:
            raise ValueError("Selected cohort labels do not cover frozen membership")
        return np.asarray([found[pair] for pair in pairs], dtype=object)

    def bind_config(self, config: dict[str, Any]) -> None:
        self.config = config

    def load_fit(self, snapshot: Snapshot) -> FitData:
        pairs = snapshot.fit_pairs
        return FitData(
            self._fit_expression(snapshot),
            self._fit_labels(snapshot),
            tuple(patient for patient, _ in pairs),
            tuple(sample for _, sample in pairs),
        )

    def load_test_expression(self, snapshot: Snapshot) -> TestData:
        self._require_access("test_access_started")
        if not self.full_provenance_verified:
            raise AuthorizationError("Test expression requires the explicit test-access stage")
        pairs = snapshot.test_pairs
        return TestData(
            self._test_expression(snapshot),
            tuple(patient for patient, _ in pairs),
            tuple(sample for _, sample in pairs),
        )

    def load_test_labels(self, snapshot: Snapshot) -> np.ndarray:
        return self._test_labels_from_cohort(snapshot.test_pairs)

    def _verify_prediction_seal(self, snapshot: Snapshot) -> None:
        approved = self._require_access("predictions_sealed")
        output = self.root / "data/processed/final_test_v1"
        seal_path = output / "prediction_seal.json"
        state_path = output / "run_state.jsonl"
        entries = [json.loads(line) for line in state_path.read_text().splitlines()]
        seal = json.loads(seal_path.read_bytes())
        prediction_sha = sha256_file(output / "test_predictions.tsv")
        seal_sha = sha256_file(seal_path)
        expected = {
            "status": "sealed_before_labels",
            "attempt_id": approved["attempt_id"],
            "authorization_sha256": approved["authorization_sha256"],
            "source_revision": approved["source_revision"],
            "model_sha256": approved["model_sha256"],
            "test_identity_sha256": identity_digest(snapshot.test_pairs),
            "prediction_sha256": prediction_sha,
            "row_count": len(snapshot.test_pairs),
            "class_order": list(LABELS),
        }
        if (
            not entries
            or entries[-1].get("state") != "predictions_sealed"
            or seal != expected
            or entries[-1].get("prediction_sha256") != prediction_sha
            or entries[-1].get("seal_sha256") != seal_sha
            or entries[-1].get("test_identity_sha256") != expected["test_identity_sha256"]
            or (self.prediction_sha256 is not None and prediction_sha != self.prediction_sha256)
            or (self.prediction_seal_sha256 is not None and seal_sha != self.prediction_seal_sha256)
        ):
            raise AuthorizationError("Prediction seal differs before label access")

    def note_predictions_saved(
        self, output: Path | None = None, snapshot: Snapshot | None = None
    ) -> None:
        if not self.full_provenance_verified:
            raise AuthorizationError("Predictions cannot be sealed before test authorization")
        if output != self.root / "data/processed/final_test_v1":
            raise AuthorizationError("Prediction seal output path differs")
        if snapshot is None:
            raise AuthorizationError("Frozen test identities are required for sealing")
        self._verify_prediction_seal(snapshot)
        seal_path = output / "prediction_seal.json"
        seal = json.loads(seal_path.read_bytes())
        self.prediction_sha256 = seal["prediction_sha256"]
        self.prediction_seal_sha256 = sha256_file(seal_path)
        self.predictions_saved = True


def verify_real_authorization(
    root: Path, config_path: Path, config: dict[str, Any], authorization_path: Path | None
) -> tuple[str, str]:
    """Gate before any real fit; a separate later approval must create the local record."""
    if authorization_path is None or not authorization_path.is_file():
        raise AuthorizationError("No separate final-test authorization record was provided")
    output = (root / config["outputs"]["directory"]).resolve()
    if output.exists():
        raise FileExistsError("Final-test output already exists; no rerun or overwrite")
    revision = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root, text=True).strip()
    dirty = subprocess.check_output(["git", "status", "--porcelain"], cwd=root, text=True).strip()
    if dirty:
        raise AuthorizationError("Real final-test run requires a clean frozen checkout")
    config_sha = sha256_file(config_path)
    record = json.loads(authorization_path.read_text(encoding="utf-8"))
    if record != {
        "authorization_schema": "final_test_single_run_v1",
        "intent": "authorize_single_final_test_execution",
        "approved_commit": revision,
        "config_sha256": config_sha,
        "split_sha256": config["inputs"]["expected_sha256"]["frozen_split"],
        "output_directory": config["outputs"]["directory"],
    }:
        raise AuthorizationError("Authorization record does not bind this exact revision and run")
    return revision, sha256_file(authorization_path)


def run_real(
    root: Path, config_path: Path, authorization_path: Path | None = None
) -> dict[str, Any]:
    """Real CLI entry point. Phase 2D-B has no authorization record or CLI test flag."""
    import yaml

    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    validate_protocol(config)
    verify_real_authorization(root, config_path, config, authorization_path)
    # A future, separately approved code revision must remove this guard.
    # Keeping it here ensures neither a synthetic fixture nor a local file can
    # accidentally authorize a real fit or test access in Phase 2D-B.
    raise AuthorizationError(
        "Real 675-patient fit and final-test access are disabled until a later approved checkpoint"
    )
