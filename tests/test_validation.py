"""Synthetic safety and reproducibility checks for the one-shot validation stage."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import yaml

from oncorna import validation
from oncorna.modeling import TrainingData, scores
from oncorna.preprocessing import sha256_file
from scripts import run_validation_v1

LABELS = ["Luminal A", "Luminal B", "Basal-like", "HER2-enriched", "Normal-like"]


def small_modeling_config() -> dict:
    return {
        "labels": LABELS,
        "expression_gt": 1.0,
        "minimum_sample_fraction": 0.20,
        "penalty": "l2",
        "solver": "lbfgs",
        "class_weight": None,
        "fit_intercept": True,
        "max_iter": 200,
        "tol": 1e-4,
    }


def small_inputs() -> validation.ValidationInputs:
    rng = np.random.default_rng(123)
    train_y = np.repeat(LABELS, 8)
    validation_y = np.repeat(LABELS, 2)
    train_values = rng.normal(size=(40, 6)) + 2
    validation_values = rng.normal(size=(10, 6)) + 2
    train_values[:, 0] = np.repeat(np.arange(5), 8) + 2
    validation_values[:, 0] = np.repeat(np.arange(5), 2) + 2
    genes = [f"gene_{i}" for i in range(6)]
    train_samples = tuple(f"train_sample_{i}" for i in range(40))
    validation_samples = tuple(f"validation_sample_{i}" for i in range(10))
    training = TrainingData(
        pd.DataFrame(train_values, index=train_samples, columns=genes),
        train_y,
        tuple(f"train_patient_{i}" for i in range(40)),
        train_samples,
        {},
    )
    return validation.ValidationInputs(
        training,
        pd.DataFrame(validation_values, index=validation_samples, columns=genes),
        validation_y,
        tuple(f"validation_patient_{i}" for i in range(10)),
        validation_samples,
        "phase2b_hash",
    )


def test_fit_uses_training_only_and_probabilities_follow_fixed_order() -> None:
    inputs, config = small_inputs(), small_modeling_config()
    pipeline, result, rows = validation.fit_and_evaluate(inputs, config, 1.0)
    assert len(rows) == len(inputs.y_validation) == 10
    assert len({row["patient_id"] for row in rows}) == 10
    assert [row["sample_id"] for row in rows] == list(inputs.validation_sample_ids)
    assert pipeline.named_steps["gene_filter"].fitting_sample_ids_ == inputs.training.sample_ids
    assert pipeline.named_steps["scaler"].n_samples_seen_ == len(inputs.training.y)
    assert result["training"]["gene_prevalence_required_count"] == 8
    assert result["training"]["selected_C"] == 1.0
    assert (
        result["validation"]["macro_f1"]
        == scores(inputs.y_validation, pipeline.predict(inputs.X_validation), LABELS)["macro_f1"]
    )
    raw_prob = pipeline.predict_proba(inputs.X_validation)
    classes = list(pipeline.named_steps["classifier"].classes_)
    for i, row in enumerate(rows):
        for label in LABELS:
            assert row[f"probability_{label}"] == pytest.approx(raw_prob[i, classes.index(label)])
        assert sum(row[f"probability_{label}"] for label in LABELS) == pytest.approx(1)
    assert result["validation"]["predicted_class_counts"] == {
        label: sum(row["predicted_label"] == label for row in rows) for label in LABELS
    }
    matrix = np.asarray(result["validation"]["confusion_matrix"])
    normalized = np.asarray(result["validation"]["row_normalized_confusion_matrix"])
    np.testing.assert_allclose(normalized, matrix / matrix.sum(axis=1, keepdims=True))


def test_held_out_expression_cannot_change_training_preprocessing() -> None:
    inputs, config = small_inputs(), small_modeling_config()
    first, first_result, first_rows = validation.fit_and_evaluate(inputs, config, 1.0)
    changed_values = inputs.X_validation.copy()
    changed_values.iloc[:, :] = 10000
    changed = validation.ValidationInputs(
        inputs.training,
        changed_values,
        inputs.y_validation,
        inputs.validation_patient_ids,
        inputs.validation_sample_ids,
        inputs.phase2b_manifest_sha256,
    )
    second, _, _ = validation.fit_and_evaluate(changed, config, 1.0)
    np.testing.assert_array_equal(
        first.named_steps["gene_filter"].get_support(),
        second.named_steps["gene_filter"].get_support(),
    )
    np.testing.assert_array_equal(
        first.named_steps["scaler"].mean_, second.named_steps["scaler"].mean_
    )
    repeated, repeated_result, repeated_rows = validation.fit_and_evaluate(inputs, config, 1.0)
    assert first_rows == repeated_rows
    assert first_result["validation"] == repeated_result["validation"]
    assert first.named_steps["classifier"].n_iter_.tolist() == (
        repeated.named_steps["classifier"].n_iter_.tolist()
    )


def test_validation_loader_indexes_only_validation_expression(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    training = small_inputs().training
    genes = list(training.X.columns)
    train_pairs = [(f"train_patient_{i}", f"train_sample_{i}") for i in range(506)]
    validation_pairs = [(f"validation_patient_{i}", f"validation_sample_{i}") for i in range(169)]
    test_pairs = [(f"test_patient_{i}", f"test_sample_{i}") for i in range(169)]
    train_samples = tuple(sample for _, sample in train_pairs)
    fake_training = TrainingData(
        pd.DataFrame(np.ones((506, 6)), index=train_samples, columns=genes),
        np.repeat(LABELS, [253, 115, 85, 40, 13]),
        tuple(patient for patient, _ in train_pairs),
        train_samples,
        {},
    )
    parts = {
        "train": tuple(train_pairs),
        "validation": tuple(validation_pairs),
        "test": tuple(test_pairs),
    }
    all_pairs = train_pairs + validation_pairs + test_pairs
    validation_labels = np.repeat(LABELS, [84, 38, 28, 14, 5])
    cohort = [
        {
            "patient_id": patient,
            "sample_id": sample,
            "pam50_normalized_label": (
                validation_labels[i - 506] if 506 <= i < 675 else "Luminal A"
            ),
        }
        for i, (patient, sample) in enumerate(all_pairs)
    ]
    (tmp_path / "split.json").write_text("{}", encoding="utf-8")
    (tmp_path / "schema.json").write_text(
        json.dumps({"gene_ids": genes, "sample_ids": [sample for _, sample in all_pairs]}),
        encoding="utf-8",
    )
    selected_indices = []

    class SpyArray:
        def __getitem__(self, key: object) -> np.ndarray:
            indices = key[0] if isinstance(key, tuple) else key
            selected_indices.extend(indices)
            return np.ones((len(indices), len(genes)))

    monkeypatch.setattr(validation, "load_training_data", lambda root, config: fake_training)
    monkeypatch.setattr(validation, "read_classification_cohort", lambda *a, **k: cohort)
    monkeypatch.setattr(validation, "validate_frozen_partitions", lambda *a: parts)
    monkeypatch.setattr(validation.np, "load", lambda *a, **k: SpyArray())
    modeling = {
        "split": "split.json",
        "cohort": "cohort.tsv",
        "all_gene_schema": "schema.json",
        "all_gene_expression": "all.npy",
        "gene_count": 6,
        "labels": LABELS,
    }
    approval = {"training_size": 506, "validation_size": 169, "test_size": 169}
    loaded = validation.load_validation_inputs(tmp_path, approval, modeling, "phase2b_hash")
    assert selected_indices == list(range(506, 675))
    assert len(loaded.y_validation) == 169
    assert set(loaded.validation_patient_ids).isdisjoint(patient for patient, _ in test_pairs)


def test_existing_output_blocks_re_evaluation_before_fitting(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    config = yaml.safe_load(
        (Path(__file__).resolve().parents[1] / "configs/validation_v1.yaml").read_text(
            encoding="utf-8"
        )
    )
    config_path = tmp_path / "approval.yaml"
    config_path.write_text(yaml.safe_dump(config), encoding="utf-8")
    (tmp_path / config["output_dir"]).mkdir(parents=True)
    monkeypatch.setattr(run_validation_v1, "ROOT", tmp_path)
    monkeypatch.setattr(
        run_validation_v1,
        "git_value",
        lambda *args: config["approved_commit"] if args[0] == "rev-parse" else "main",
    )
    with pytest.raises(ValueError, match="new local-only versioned directory"):
        run_validation_v1.preflight(config_path)


def test_phase2b_result_hash_and_selection_are_verified(tmp_path: Path) -> None:
    modeling_path = tmp_path / "modeling.yaml"
    modeling_path.write_text("locked: true\n", encoding="utf-8")
    result_path = tmp_path / "candidate_results.json"
    result_path.write_text(
        json.dumps([{"C": 0.1, "fold_macro_f1_mean": 0.4}, {"C": 1.0, "fold_macro_f1_mean": 0.5}]),
        encoding="utf-8",
    )
    manifest_path = tmp_path / "manifest.json"
    manifest_path.write_text(
        json.dumps(
            {
                "version": "modeling_v1",
                "scope": "training_only_cv",
                "config_sha256": sha256_file(modeling_path),
                "selected_C": 1.0,
                "training_count": 506,
                "class_order": LABELS,
                "output_hashes": {"candidate_results.json": sha256_file(result_path)},
            }
        ),
        encoding="utf-8",
    )
    approval = {
        "modeling_config": "modeling.yaml",
        "modeling_config_sha256": sha256_file(modeling_path),
        "phase2b_manifest": "manifest.json",
        "phase2b_results": "candidate_results.json",
        "selected_C": 1.0,
        "training_size": 506,
    }
    assert validation.verify_phase2b(tmp_path, approval, {"labels": LABELS}) == sha256_file(
        manifest_path
    )
    result_path.write_text("[]", encoding="utf-8")
    with pytest.raises(ValueError, match="checksum"):
        validation.verify_phase2b(tmp_path, approval, {"labels": LABELS})
