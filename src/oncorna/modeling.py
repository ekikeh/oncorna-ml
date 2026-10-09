"""Training-only, patient-stratified PAM50 cross-validation."""

from __future__ import annotations

import json
import time
import warnings
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from sklearn.dummy import DummyClassifier
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    accuracy_score,
    balanced_accuracy_score,
    confusion_matrix,
    precision_recall_fscore_support,
)
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from oncorna.preprocessing import (
    GenePrevalenceFilter,
    sha256_file,
    validate_frozen_partitions,
    verify_split_provenance,
)
from oncorna.qc import read_classification_cohort


@dataclass(frozen=True)
class TrainingData:
    X: pd.DataFrame
    y: np.ndarray
    patient_ids: tuple[str, ...]
    sample_ids: tuple[str, ...]
    input_hashes: dict[str, str]


def _verified_path(root: Path, relative: str, expected: str) -> Path:
    path = root / relative
    actual = sha256_file(path)
    if actual != expected:
        raise ValueError(f"Input checksum mismatch: {relative}")
    return path


def load_training_data(root: Path, config: dict[str, Any]) -> TrainingData:
    """Verify frozen inputs, then materialize only training rows of expression."""
    paths = {
        name: _verified_path(root, config[name], config[f"expected_{suffix}_sha256"])
        for name, suffix in (
            ("split", "split"),
            ("cohort", "cohort"),
            ("source_expression", "source"),
            ("all_gene_expression", "all_gene"),
            ("all_gene_schema", "schema"),
        )
    }
    manifest_path = root / config["preprocessing_manifest"]
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("version") != "preprocessing_v1":
        raise ValueError("Unexpected preprocessing manifest version")
    for section, key, expected in (
        ("split", "sha256", config["expected_split_sha256"]),
        ("cohort", "sha256", config["expected_cohort_sha256"]),
        ("source", "sha256", config["expected_source_sha256"]),
        ("all_gene_view", "sha256", config["expected_all_gene_sha256"]),
        ("all_gene_view", "schema_sha256", config["expected_schema_sha256"]),
    ):
        if manifest[section][key] != expected:
            raise ValueError(f"Preprocessing manifest disagrees with {section}.{key}")
    split = json.loads(paths["split"].read_text(encoding="utf-8"))
    cohort = read_classification_cohort(paths["cohort"], expected_sample_count=844)
    partitions = validate_frozen_partitions(cohort, split)
    verify_split_provenance(
        split,
        cohort_path=paths["cohort"],
        historical_config_path=root / "configs/default.yaml",
        historical_matrix_path=root / "data/processed/ml_expression_matrix.tsv.gz",
        historical_matrix_manifest_path=root / "data/processed/ml_matrix_manifest.json",
    )
    train_pairs = partitions["train"]
    if len(train_pairs) != config["training_size"]:
        raise ValueError("Training size differs from approved configuration")
    train_patients = {patient for patient, _ in train_pairs}
    held_out_patients = {
        patient for name in ("validation", "test") for patient, _ in partitions[name]
    }
    if train_patients & held_out_patients:
        raise ValueError("Held-out patient entered training partition")
    schema = json.loads(paths["all_gene_schema"].read_text(encoding="utf-8"))
    sample_ids = tuple(schema["sample_ids"])
    gene_ids = tuple(schema["gene_ids"])
    if (
        schema["version"] != "preprocessing_v1"
        or schema["orientation"] != "samples_by_genes"
        or schema["expression_units"] != "log2(normalized_count + 1)"
        or schema["dtype"] != "float64"
        or schema["shape"] != [844, config["gene_count"]]
        or len(set(gene_ids)) != len(gene_ids)
        or len(gene_ids) != config["gene_count"]
        or sample_ids != tuple(row["sample_id"] for row in cohort)
    ):
        raise ValueError("All-gene schema disagrees with approved cohort or configuration")
    values = np.load(paths["all_gene_expression"], mmap_mode="r", allow_pickle=False)
    if values.shape != tuple(schema["shape"]) or values.dtype != np.float64:
        raise ValueError("All-gene array disagrees with schema")
    positions = {sample: index for index, sample in enumerate(sample_ids)}
    training_samples = tuple(sample for _, sample in train_pairs)
    selected = values[[positions[sample] for sample in training_samples], :]
    if not np.isfinite(selected).all():
        raise ValueError("Training expression has nonfinite values; no imputer approved")
    X = pd.DataFrame(selected, index=training_samples, columns=gene_ids)
    labels_by_pair = {
        (row["patient_id"], row["sample_id"]): row["pam50_normalized_label"] for row in cohort
    }
    y = np.array([labels_by_pair[pair] for pair in train_pairs], dtype=object)
    if Counter(y) != Counter(
        {
            "Luminal A": 253,
            "Luminal B": 115,
            "Basal-like": 85,
            "HER2-enriched": 40,
            "Normal-like": 13,
        }
    ):
        raise ValueError("Training class counts differ from frozen split")
    if set(y) != set(config["labels"]):
        raise ValueError("Class labels differ from approved configuration")
    return TrainingData(
        X=X,
        y=y,
        patient_ids=tuple(patient for patient, _ in train_pairs),
        sample_ids=training_samples,
        input_hashes={
            **{name: sha256_file(path) for name, path in paths.items()},
            "preprocessing_manifest": sha256_file(manifest_path),
        },
    )


def make_folds(data: TrainingData, n_splits: int, seed: int) -> list[dict[str, Any]]:
    if len(set(data.patient_ids)) != len(data.patient_ids):
        raise ValueError("Patient IDs are not unique")
    splitter = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    folds = []
    held_out_indices = []
    for number, (fit, hold) in enumerate(splitter.split(data.X, data.y), start=1):
        fit_patients = [data.patient_ids[index] for index in fit]
        hold_patients = [data.patient_ids[index] for index in hold]
        if set(fit_patients) & set(hold_patients):
            raise ValueError("Patient overlaps within CV fold")
        held_out_indices.extend(hold.tolist())
        folds.append(
            {
                "fold": number,
                "fit_indices": fit.tolist(),
                "hold_indices": hold.tolist(),
                "fit_patient_ids": fit_patients,
                "hold_patient_ids": hold_patients,
                "fit_class_counts": dict(Counter(data.y[fit])),
                "hold_class_counts": dict(Counter(data.y[hold])),
            }
        )
    if sorted(held_out_indices) != list(range(len(data.patient_ids))):
        raise ValueError("Every training patient must occur in one held-out fold")
    return folds


def make_pipeline(config: dict[str, Any], c_value: float) -> Pipeline:
    return Pipeline(
        [
            (
                "gene_filter",
                GenePrevalenceFilter(config["expression_gt"], config["minimum_sample_fraction"]),
            ),
            ("scaler", StandardScaler()),
            (
                "classifier",
                LogisticRegression(
                    C=c_value,
                    penalty=config["penalty"],
                    solver=config["solver"],
                    class_weight=config["class_weight"],
                    fit_intercept=config["fit_intercept"],
                    max_iter=config["max_iter"],
                    tol=config["tol"],
                ),
            ),
        ]
    )


def scores(y_true: np.ndarray, y_pred: np.ndarray, labels: list[str]) -> dict[str, Any]:
    precision, recall, f1, support = precision_recall_fscore_support(
        y_true,
        y_pred,
        labels=labels,
        zero_division=0,
    )
    return {
        "macro_f1": float(np.mean(f1)),
        "balanced_accuracy": float(balanced_accuracy_score(y_true, y_pred)),
        "accuracy": float(accuracy_score(y_true, y_pred)),
        "per_class": {
            label: {
                "precision": float(precision[i]),
                "recall": float(recall[i]),
                "f1": float(f1[i]),
                "support": int(support[i]),
            }
            for i, label in enumerate(labels)
        },
        "confusion_matrix": confusion_matrix(y_true, y_pred, labels=labels).tolist(),
    }


def evaluate_candidate(
    data: TrainingData,
    folds: list[dict[str, Any]],
    config: dict[str, Any],
    c_value: float | None,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    """Run fresh estimators per fold and return aggregate and local-only predictions."""
    label_order = config["labels"]
    oof = np.empty(len(data.y), dtype=object)
    oof[:] = None
    fold_results = []
    for fold in folds:
        fit, hold = fold["fit_indices"], fold["hold_indices"]
        estimator = (
            DummyClassifier(strategy="most_frequent")
            if c_value is None
            else make_pipeline(config, c_value)
        )
        started = time.perf_counter()
        try:
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                estimator.fit(data.X.iloc[fit], data.y[fit])
        except Exception as exc:
            model = "dummy" if c_value is None else f"logistic C={c_value:g}"
            raise RuntimeError(f"{model} failed in fold {fold['fold']}: {exc}") from exc
        elapsed = time.perf_counter() - started
        pred = estimator.predict(data.X.iloc[hold])
        oof[hold] = pred
        warning_records = [
            {"category": item.category.__name__, "message": str(item.message)} for item in caught
        ]
        convergence = any(issubclass(item.category, ConvergenceWarning) for item in caught)
        gene_count = None
        iterations = None
        if c_value is not None:
            gene_filter = estimator.named_steps["gene_filter"]
            scaler = estimator.named_steps["scaler"]
            if gene_filter.fitting_sample_ids_ != tuple(data.sample_ids[i] for i in fit):
                raise ValueError("Gene filter fit rows differ from CV fitting rows")
            if np.any(np.asarray(scaler.n_samples_seen_) != len(fit)):
                raise ValueError("Scaler fit count differs from CV fitting count")
            gene_count = len(gene_filter.selected_gene_ids_)
            iterations = estimator.named_steps["classifier"].n_iter_.tolist()
        fold_results.append(
            {
                "fold": fold["fold"],
                "fit_size": len(fit),
                "hold_size": len(hold),
                "retained_gene_count": gene_count,
                "fit_seconds": elapsed,
                "convergence_warning": convergence,
                "warnings": warning_records,
                "classifier_iterations": iterations,
                "metrics": scores(data.y[hold], pred, label_order),
            }
        )
    if any(value is None for value in oof):
        raise ValueError("Missing out-of-fold prediction")
    macros = [item["metrics"]["macro_f1"] for item in fold_results]
    result = {
        "model": "dummy_most_frequent" if c_value is None else "logistic_l2_multinomial",
        "C": c_value,
        "folds": fold_results,
        "fold_macro_f1_mean": float(np.mean(macros)),
        "fold_macro_f1_std": float(np.std(macros, ddof=1)),
        "fold_macro_f1_range": [float(min(macros)), float(max(macros))],
        "pooled_oof_metrics": scores(data.y, oof, label_order),
        "any_convergence_warning": any(item["convergence_warning"] for item in fold_results),
        "total_fit_seconds": float(sum(item["fit_seconds"] for item in fold_results)),
    }
    predictions = [
        {
            "patient_id": data.patient_ids[i],
            "sample_id": data.sample_ids[i],
            "true_label": str(data.y[i]),
            "predicted_label": str(oof[i]),
            "model": result["model"],
            "C": c_value,
        }
        for i in range(len(data.y))
    ]
    return result, predictions


def select_c(candidate_results: list[dict[str, Any]]) -> float:
    candidates = [result for result in candidate_results if result["C"] is not None]
    if not candidates or any(not np.isfinite(item["fold_macro_f1_mean"]) for item in candidates):
        raise ValueError("No valid logistic candidate scores")
    return float(
        sorted(candidates, key=lambda item: (-item["fold_macro_f1_mean"], item["C"]))[0]["C"]
    )
