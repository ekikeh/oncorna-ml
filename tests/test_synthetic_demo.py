"""Standalone demo regression tests; all values and identities are invented."""

from __future__ import annotations

import builtins
import json
import pickle
import socket
import subprocess
import sys
from dataclasses import replace
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
import pytest
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

import oncorna.synthetic_demo as demo
from oncorna.prediction import PredictionInputError


def test_generator_is_repeatable_and_identifiers_are_invented() -> None:
    first = demo.generate_synthetic_dataset()
    second = demo.generate_synthetic_dataset()
    for attribute in ("train_X", "evaluation_X"):
        pd.testing.assert_frame_equal(getattr(first, attribute), getattr(second, attribute))
    for attribute in ("train_y", "evaluation_y"):
        pd.testing.assert_series_equal(getattr(first, attribute), getattr(second, attribute))
    assert first.seed == demo.DEMO_SEED == 20261026
    assert first.train_X.shape == (100, 12)
    assert first.evaluation_X.shape == (40, 12)
    assert tuple(first.train_X.columns) == tuple(first.evaluation_X.columns) == demo.DEMO_GENES
    assert set(first.train_X.index).isdisjoint(first.evaluation_X.index)
    assert set(first.train_X.index) | set(first.evaluation_X.index) == {
        f"SYN-BRCA-{index:04d}" for index in range(1, 141)
    }
    for expression, labels, support in (
        (first.train_X, first.train_y, 20),
        (first.evaluation_X, first.evaluation_y, 8),
    ):
        assert np.isfinite(expression.to_numpy()).all()
        assert expression.to_numpy().min() >= 0
        assert expression.index.equals(labels.index)
        assert labels.value_counts().to_dict() == dict.fromkeys(demo.DEMO_CLASSES, support)
    assert not first.train_X.equals(demo.generate_synthetic_dataset(7).train_X)


@pytest.mark.parametrize("seed", [True, -1, 2**32, 1.5, "1", None])
def test_invalid_seed_is_rejected(seed: object) -> None:
    with pytest.raises(ValueError, match="Seed"):
        demo.generate_synthetic_dataset(seed)


def test_repeated_predictions_and_cli_results_are_identical(capsys: pytest.CaptureFixture) -> None:
    results = []
    original_predict = Pipeline.predict

    def record_predictions(self: Pipeline, *args: object, **kwargs: object):
        predictions = original_predict(self, *args, **kwargs)
        results.append(predictions.copy())
        return predictions

    from unittest.mock import patch

    with patch.object(Pipeline, "predict", record_predictions):
        first = demo.run_demo()
        second = demo.run_demo()
    np.testing.assert_array_equal(results[0], results[1])
    assert first == second
    assert first["scope"] == "synthetic_demo_only"
    assert first["model"]["source"] == "fresh_synthetic_training_only"
    assert first["training_samples"] == 100 and first["evaluation_samples"] == 40
    assert first["class_order"] == list(demo.DEMO_CLASSES)
    assert first["confusion_matrix"]["rows"] == "true_synthetic_class"
    assert first["confusion_matrix"]["columns"] == "predicted_synthetic_class"
    matrix = np.asarray(first["confusion_matrix"]["counts"])
    assert matrix.shape == (5, 5) and matrix.sum() == 40
    np.testing.assert_array_equal(matrix.sum(axis=1), np.repeat(8, 5))
    assert np.isfinite(first["metrics"]["macro_f1"])
    assert first["training_data_sha256"] != first["evaluation_data_sha256"]
    assert demo.main([]) == 0
    assert json.loads(capsys.readouterr().out) == first


def test_only_training_fits_preprocessing_and_prediction_never_refits(monkeypatch) -> None:
    data = demo.generate_synthetic_dataset()
    fit_calls = []
    predictions = []
    original_fit, original_predict = Pipeline.fit, Pipeline.predict

    def observe_fit(self: Pipeline, expression: pd.DataFrame, labels: pd.Series, **kwargs):
        fit_calls.append(tuple(expression.index))
        return original_fit(self, expression, labels, **kwargs)

    def observe_prediction(self: Pipeline, expression: pd.DataFrame, **kwargs):
        scaler = self.named_steps["scaler"]
        assert scaler.n_samples_seen_ == 100
        np.testing.assert_allclose(scaler.mean_, data.train_X.mean().to_numpy())
        means, coefficients = scaler.mean_.copy(), self.named_steps["classifier"].coef_.copy()

        def forbidden_fit(*_args, **_kwargs):
            pytest.fail("Fitting was attempted during prediction")

        with monkeypatch.context() as guard:
            guard.setattr(Pipeline, "fit", forbidden_fit)
            guard.setattr(StandardScaler, "fit", forbidden_fit)
            result = original_predict(self, expression, **kwargs)
        np.testing.assert_array_equal(scaler.mean_, means)
        np.testing.assert_array_equal(self.named_steps["classifier"].coef_, coefficients)
        predictions.append(result)
        return result

    monkeypatch.setattr(Pipeline, "fit", observe_fit)
    monkeypatch.setattr(Pipeline, "predict", observe_prediction)
    demo.run_demo()
    assert fit_calls == [tuple(data.train_X.index)]
    assert len(predictions) == 1


def test_evaluation_values_do_not_change_fitted_state(monkeypatch) -> None:
    data = demo.generate_synthetic_dataset()
    states = []
    original_fit = Pipeline.fit

    def observe_fit(self: Pipeline, *args, **kwargs):
        result = original_fit(self, *args, **kwargs)
        states.append(
            (self.named_steps["scaler"].mean_.copy(), self.named_steps["classifier"].coef_.copy())
        )
        return result

    monkeypatch.setattr(Pipeline, "fit", observe_fit)
    demo.run_demo()
    changed = data.evaluation_X * 100.0
    monkeypatch.setattr(
        demo, "generate_synthetic_dataset", lambda _seed: replace(data, evaluation_X=changed)
    )
    demo.run_demo()
    for first, second in zip(states[0], states[1], strict=True):
        np.testing.assert_array_equal(first, second)


def test_cli_validation_error_has_nonzero_exit_and_clear_message(monkeypatch, capsys) -> None:
    data = demo.generate_synthetic_dataset()
    malformed = replace(data, evaluation_X=data.evaluation_X.iloc[:, :-1])
    monkeypatch.setattr(demo, "generate_synthetic_dataset", lambda _seed: malformed)
    with pytest.raises(PredictionInputError, match="missing, extra or reordered"):
        demo.run_demo()
    with pytest.raises(SystemExit) as error:
        demo.main([])
    assert error.value.code == 2
    output = capsys.readouterr()
    assert not output.out
    assert "Synthetic demo input validation failed" in output.err


@pytest.mark.parametrize("option", ["--model", "--input", "--output", "--authorization"])
def test_cli_rejects_artifact_and_path_options(option: str, capsys) -> None:
    with pytest.raises(SystemExit) as error:
        demo.main([option, "SYN-NOT-A-REAL-PATH"])
    assert error.value.code == 2
    assert "unrecognized arguments" in capsys.readouterr().err


def test_no_caller_model_or_fallback_can_be_substituted(monkeypatch) -> None:
    with pytest.raises(TypeError, match="unexpected keyword"):
        demo.run_demo(model=object())

    def fail_fit(*_args, **_kwargs):
        raise RuntimeError("Synthetic fitting failed; no fallback allowed")

    monkeypatch.setattr(Pipeline, "fit", fail_fit)
    with pytest.raises(RuntimeError, match="no fallback"):
        demo.run_demo()


def test_demo_never_reads_files_loads_models_or_uses_network(monkeypatch) -> None:
    # Imported dependencies are warm; any runtime filesystem/network/model access fails.
    def forbidden(*_args, **_kwargs):
        pytest.fail("Standalone synthetic demo attempted file, model or network access")

    with monkeypatch.context() as guard:
        guard.setattr(builtins, "open", forbidden)
        guard.setattr(Path, "open", forbidden)
        guard.setattr(Path, "read_bytes", forbidden)
        guard.setattr(Path, "read_text", forbidden)
        guard.setattr(joblib, "load", forbidden)
        guard.setattr(pickle, "load", forbidden)
        guard.setattr(np, "load", forbidden)
        guard.setattr(socket, "socket", forbidden)
        guard.setattr(subprocess, "Popen", forbidden)
        result = demo.run_demo()
    assert result["scope"] == demo.SCOPE


def test_fresh_cli_process_does_not_import_protected_evaluators() -> None:
    code = (
        "import sys; from oncorna.synthetic_demo import run_demo; "
        "assert run_demo()['scope'] == 'synthetic_demo_only'; "
        "assert not any(name in sys.modules for name in "
        "['oncorna.final_test', 'oncorna.final_test_consumer', "
        "'oncorna.final_model_fit', 'oncorna.final_model_fit_real', "
        "'oncorna.final_fit_real', 'oncorna.windows_attempt_claim', 'oncorna.synthetic_paths'])"
    )
    subprocess.run([sys.executable, "-B", "-c", code], check=True, capture_output=True, text=True)
