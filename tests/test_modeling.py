"""Synthetic checks for the Phase 2B training-only CV boundary."""

from __future__ import annotations

import json
from copy import deepcopy
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import yaml
from sklearn.exceptions import ConvergenceWarning

from oncorna.modeling import (
    TrainingData,
    evaluate_candidate,
    make_folds,
    make_pipeline,
    scores,
    select_c,
)
from scripts import run_training_cv

LABELS = ["Luminal A", "Luminal B", "Basal-like", "HER2-enriched", "Normal-like"]


def fixture_data() -> TrainingData:
    labels = np.repeat(LABELS, 8)
    generator = np.random.default_rng(91)
    values = generator.normal(size=(40, 6))
    values[:, 0] = np.repeat(np.arange(5), 8) + generator.normal(0, 0.1, 40)
    values[:, 1] = 2.0
    samples = tuple(f"training_sample_{i}" for i in range(40))
    patients = tuple(f"training_patient_{i}" for i in range(40))
    frame = pd.DataFrame(values, index=samples, columns=[f"gene_{i}" for i in range(6)])
    return TrainingData(frame, labels, patients, samples, {})


def fixture_config() -> dict:
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


def test_folds_are_patient_disjoint_and_cover_training_once() -> None:
    data = fixture_data()
    folds = make_folds(data, 4, 20261009)
    assert folds == make_folds(data, 4, 20261009)
    assert sorted(index for fold in folds for index in fold["hold_indices"]) == list(range(40))
    for fold in folds:
        assert not set(fold["fit_patient_ids"]) & set(fold["hold_patient_ids"])
        assert set(fold["fit_patient_ids"] + fold["hold_patient_ids"]) == set(data.patient_ids)
        assert fold["hold_class_counts"] == dict.fromkeys(LABELS, 2)
    duplicate = TrainingData(data.X, data.y, (data.patient_ids[0],) * 40, data.sample_ids, {})
    with pytest.raises(ValueError, match="not unique"):
        make_folds(duplicate, 4, 20261009)


def test_oof_predictions_score_order_and_repeated_execution() -> None:
    data, config = fixture_data(), fixture_config()
    folds = make_folds(data, 4, 20261009)
    first, predictions = evaluate_candidate(data, folds, config, 0.01)
    second, again = evaluate_candidate(data, folds, config, 0.01)
    assert len(predictions) == len(again) == len(data.y)
    assert predictions == again
    assert [row["patient_id"] for row in predictions] == list(data.patient_ids)
    assert all(row["model"] == "logistic_l2_multinomial" for row in predictions)
    assert [fold["metrics"] for fold in first["folds"]] == [
        fold["metrics"] for fold in second["folds"]
    ]
    assert list(first["pooled_oof_metrics"]["per_class"]) == LABELS
    assert sum(first["pooled_oof_metrics"]["per_class"][label]["support"] for label in LABELS) == 40
    assert all(fold["retained_gene_count"] is not None for fold in first["folds"])
    dummy, dummy_predictions = evaluate_candidate(data, folds, config, None)
    assert len(dummy_predictions) == 40
    assert dummy["model"] == "dummy_most_frequent"


def test_held_out_expression_does_not_change_fitted_preprocessing() -> None:
    data, config = fixture_data(), fixture_config()
    fold = make_folds(data, 4, 20261009)[0]
    fit, hold = fold["fit_indices"], fold["hold_indices"]
    first = make_pipeline(config, 0.01).fit(data.X.iloc[fit], data.y[fit])
    changed = data.X.copy()
    changed.iloc[hold, :] = 9999
    second = make_pipeline(config, 0.01).fit(changed.iloc[fit], data.y[fit])
    first_filter = first.named_steps["gene_filter"]
    second_filter = second.named_steps["gene_filter"]
    assert first_filter.fitting_sample_ids_ == tuple(data.sample_ids[i] for i in fit)
    assert first_filter.fitting_sample_ids_ == second_filter.fitting_sample_ids_
    np.testing.assert_array_equal(first_filter.get_support(), second_filter.get_support())
    np.testing.assert_array_equal(
        first.named_steps["scaler"].mean_, second.named_steps["scaler"].mean_
    )
    assert first.named_steps["scaler"].n_samples_seen_ == len(fit)


def test_fixed_label_mapping_and_c_tie_break() -> None:
    metric = scores(
        np.array(["Normal-like", "Luminal A"]), np.array(["Luminal A", "Luminal A"]), LABELS
    )
    assert list(metric["per_class"]) == LABELS
    assert metric["confusion_matrix"][4][0] == 1
    assert metric["per_class"]["Normal-like"]["recall"] == 0
    candidates = [
        {"C": 1.0, "fold_macro_f1_mean": 0.6},
        {"C": 0.01, "fold_macro_f1_mean": 0.6},
        {"C": 0.001, "fold_macro_f1_mean": 0.5},
    ]
    assert select_c(candidates) == 0.01


def test_convergence_warning_is_recorded() -> None:
    data, config = fixture_data(), fixture_config()
    config = deepcopy(config)
    config["max_iter"] = 1
    result, _ = evaluate_candidate(data, make_folds(data, 4, 20261009), config, 1.0)
    assert result["any_convergence_warning"]
    assert any(
        warning["category"] == ConvergenceWarning.__name__
        for fold in result["folds"]
        for warning in fold["warnings"]
    )


def test_failed_fit_writes_failure_record(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    source = Path(__file__).resolve().parents[1] / "configs/modeling_v1.yaml"
    config = yaml.safe_load(source.read_text(encoding="utf-8"))
    config_path = tmp_path / "modeling_v1.yaml"
    config_path.write_text(yaml.safe_dump(config), encoding="utf-8")
    monkeypatch.setattr(run_training_cv, "ROOT", tmp_path)
    monkeypatch.setattr(run_training_cv, "load_training_data", lambda root, cfg: fixture_data())

    def fail_candidate(*args: object) -> None:
        raise RuntimeError("synthetic fit failure")

    monkeypatch.setattr(run_training_cv, "evaluate_candidate", fail_candidate)
    with pytest.raises(RuntimeError, match="synthetic fit failure"):
        run_training_cv.run(config_path)
    failure = json.loads((tmp_path / config["output_dir"] / "failure.json").read_text())
    assert failure["candidate"] == "dummy"
    assert failure["error_type"] == "RuntimeError"
    assert failure["completed_candidates"] == []
