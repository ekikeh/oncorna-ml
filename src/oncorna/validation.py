"""Frozen, single-use validation evaluation for the approved Phase 2B model."""

from __future__ import annotations

import hashlib
import json
import math
import time
import warnings
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from sklearn.exceptions import ConvergenceWarning

from oncorna.modeling import TrainingData, load_training_data, make_pipeline, scores, select_c
from oncorna.preprocessing import sha256_file, validate_frozen_partitions
from oncorna.qc import read_classification_cohort


@dataclass(frozen=True)
class ValidationInputs:
    training: TrainingData
    X_validation: pd.DataFrame
    y_validation: np.ndarray
    validation_patient_ids: tuple[str, ...]
    validation_sample_ids: tuple[str, ...]
    phase2b_manifest_sha256: str


def verify_phase2b(root: Path, approval: dict[str, Any], modeling: dict[str, Any]) -> str:
    """Pin model selection to the completed, checksummed Phase 2B run."""
    config_path = root / approval["modeling_config"]
    if sha256_file(config_path) != approval["modeling_config_sha256"]:
        raise ValueError("Phase 2B modeling configuration changed")
    manifest_path = root / approval["phase2b_manifest"]
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if (
        manifest["version"] != "modeling_v1"
        or manifest["scope"] != "training_only_cv"
        or manifest["config_sha256"] != approval["modeling_config_sha256"]
        or manifest["selected_C"] != approval["selected_C"]
        or manifest["training_count"] != approval["training_size"]
        or manifest["class_order"] != modeling["labels"]
    ):
        raise ValueError("Phase 2B manifest disagrees with approved selection")
    output_dir = manifest_path.parent
    for name, expected in manifest["output_hashes"].items():
        if sha256_file(output_dir / name) != expected:
            raise ValueError(f"Phase 2B output checksum mismatch: {name}")
    if (
        sha256_file(root / approval["phase2b_results"])
        != manifest["output_hashes"]["candidate_results.json"]
    ):
        raise ValueError("Phase 2B result path differs from manifest")
    results = json.loads((root / approval["phase2b_results"]).read_text(encoding="utf-8"))
    if select_c(results) != approval["selected_C"]:
        raise ValueError("Phase 2B candidate scores disagree with selected C")
    return sha256_file(manifest_path)


def load_validation_inputs(
    root: Path,
    approval: dict[str, Any],
    modeling: dict[str, Any],
    phase2b_manifest_sha256: str,
) -> ValidationInputs:
    """Verify all partition identities; materialize training and validation rows only."""
    training = load_training_data(root, modeling)
    split = json.loads((root / modeling["split"]).read_text(encoding="utf-8"))
    cohort = read_classification_cohort(root / modeling["cohort"], expected_sample_count=844)
    parts = validate_frozen_partitions(cohort, split)
    if {name: len(parts[name]) for name in parts} != {
        "train": approval["training_size"],
        "validation": approval["validation_size"],
        "test": approval["test_size"],
    }:
        raise ValueError("Frozen partition sizes differ from Phase 2C approval")
    train_ids = set(training.patient_ids)
    validation_pairs = parts["validation"]
    validation_ids = {patient for patient, _ in validation_pairs}
    test_ids = {patient for patient, _ in parts["test"]}
    if train_ids & validation_ids or train_ids & test_ids or validation_ids & test_ids:
        raise ValueError("Frozen patient partitions overlap")
    if tuple(zip(training.patient_ids, training.sample_ids, strict=True)) != parts["train"]:
        raise ValueError("Training rows differ from frozen patient/sample pairs")
    schema = json.loads((root / modeling["all_gene_schema"]).read_text(encoding="utf-8"))
    if tuple(schema["gene_ids"]) != tuple(training.X.columns):
        raise ValueError("Validation gene schema differs from training feature order")
    sample_positions = {sample: position for position, sample in enumerate(schema["sample_ids"])}
    validation_samples = tuple(sample for _, sample in validation_pairs)
    values = np.load(root / modeling["all_gene_expression"], mmap_mode="r", allow_pickle=False)
    selected = values[[sample_positions[sample] for sample in validation_samples], :]
    if selected.shape != (approval["validation_size"], modeling["gene_count"]):
        raise ValueError("Validation expression shape differs from approved design")
    if not np.isfinite(selected).all():
        raise ValueError("Validation expression has nonfinite values; no imputer approved")
    X_validation = pd.DataFrame(selected, index=validation_samples, columns=training.X.columns)
    labels_by_pair = {
        (row["patient_id"], row["sample_id"]): row["pam50_normalized_label"] for row in cohort
    }
    y_validation = np.asarray([labels_by_pair[pair] for pair in validation_pairs], dtype=object)
    if Counter(y_validation) != Counter(
        {"Luminal A": 84, "Luminal B": 38, "Basal-like": 28, "HER2-enriched": 14, "Normal-like": 5}
    ):
        raise ValueError("Validation class counts differ from frozen split")
    if set(y_validation) != set(modeling["labels"]):
        raise ValueError("Validation labels differ from approved class order")
    return ValidationInputs(
        training=training,
        X_validation=X_validation,
        y_validation=y_validation,
        validation_patient_ids=tuple(patient for patient, _ in validation_pairs),
        validation_sample_ids=validation_samples,
        phase2b_manifest_sha256=phase2b_manifest_sha256,
    )


def ordered_gene_digest(gene_ids: tuple[str, ...]) -> str:
    payload = json.dumps(gene_ids, ensure_ascii=False, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def fit_and_evaluate(
    inputs: ValidationInputs, modeling: dict[str, Any], selected_c: float
) -> tuple[Any, dict[str, Any], list[dict[str, Any]]]:
    """Fit once on training only, then evaluate the frozen validation rows once."""
    if set(inputs.training.patient_ids) & set(inputs.validation_patient_ids):
        raise ValueError("Validation patient entered training")
    if tuple(inputs.training.X.columns) != tuple(inputs.X_validation.columns):
        raise ValueError("Training/validation feature schemas differ")
    pipeline = make_pipeline(modeling, selected_c)
    started = time.perf_counter()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        pipeline.fit(inputs.training.X, inputs.training.y)
    fit_seconds = time.perf_counter() - started
    gene_filter = pipeline.named_steps["gene_filter"]
    scaler = pipeline.named_steps["scaler"]
    classifier = pipeline.named_steps["classifier"]
    if gene_filter.fitting_sample_ids_ != inputs.training.sample_ids:
        raise ValueError("Gene filter was not fitted on exactly the frozen training samples")
    if np.any(np.asarray(scaler.n_samples_seen_) != len(inputs.training.y)):
        raise ValueError("Scaler fitting count differs from frozen training count")
    expected_required_count = math.ceil(
        modeling["minimum_sample_fraction"] * len(inputs.training.y)
    )
    if gene_filter.required_count_ != expected_required_count:
        raise ValueError("Training gene prevalence threshold count differs from approved rule")
    if tuple(classifier.classes_) != tuple(sorted(modeling["labels"])):
        raise ValueError("Classifier classes differ from approved class set")
    predicted = pipeline.predict(inputs.X_validation)
    raw_probabilities = pipeline.predict_proba(inputs.X_validation)
    class_positions = [list(classifier.classes_).index(label) for label in modeling["labels"]]
    probabilities = raw_probabilities[:, class_positions]
    if (
        len(predicted) != len(inputs.y_validation)
        or probabilities.shape != (len(inputs.y_validation), len(modeling["labels"]))
        or not np.isfinite(probabilities).all()
        or np.any(probabilities < 0)
        or np.any(probabilities > 1)
        or not np.allclose(probabilities.sum(axis=1), 1, atol=1e-9)
    ):
        raise ValueError("Validation predictions or aligned probabilities are invalid")
    metric = scores(inputs.y_validation, predicted, modeling["labels"])
    matrix = np.asarray(metric["confusion_matrix"])
    row_totals = matrix.sum(axis=1, keepdims=True)
    metric["row_normalized_confusion_matrix"] = (matrix / row_totals).tolist()
    metric["predicted_class_counts"] = {
        label: int(np.count_nonzero(predicted == label)) for label in modeling["labels"]
    }
    warnings_list = [
        {"category": item.category.__name__, "message": str(item.message)} for item in caught
    ]
    training = {
        "training_patient_count": len(inputs.training.y),
        "fit_sample_ids_sha256": ordered_gene_digest(inputs.training.sample_ids),
        "retained_gene_count": len(gene_filter.selected_gene_ids_),
        "ordered_gene_list_sha256": ordered_gene_digest(gene_filter.selected_gene_ids_),
        "gene_prevalence_required_count": gene_filter.required_count_,
        "scaler_fit_count": int(np.min(np.asarray(scaler.n_samples_seen_))),
        "solver_iterations": classifier.n_iter_.tolist(),
        "convergence_warning": any(
            issubclass(item.category, ConvergenceWarning) for item in caught
        ),
        "warnings": warnings_list,
        "fit_seconds": fit_seconds,
        "selected_C": selected_c,
    }
    rows = [
        {
            "patient_id": inputs.validation_patient_ids[i],
            "sample_id": inputs.validation_sample_ids[i],
            "true_label": str(inputs.y_validation[i]),
            "predicted_label": str(predicted[i]),
            **{
                f"probability_{label}": float(probabilities[i, j])
                for j, label in enumerate(modeling["labels"])
            },
        }
        for i in range(len(predicted))
    ]
    return pipeline, {"training": training, "validation": metric}, rows
