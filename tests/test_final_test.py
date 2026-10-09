"""Synthetic-only checks for the frozen final-test infrastructure."""

from __future__ import annotations

import copy
import csv
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import yaml

from oncorna.final_test import (
    AuthorizationError,
    FitData,
    RealSource,
    SyntheticFixtureSource,
    _score,
    execute_once,
    run_real,
    validate_memberships,
    validate_protocol,
)
from oncorna.final_test import (
    TestData as HeldOutData,
)

LABELS = ("Luminal A", "Luminal B", "Basal-like", "HER2-enriched", "Normal-like")


def synthetic_config() -> dict:
    source = Path(__file__).resolve().parents[1] / "configs/final_test_v1.yaml"
    config = copy.deepcopy(yaml.safe_load(source.read_text(encoding="utf-8")))
    config["partitions"].update(
        {
            "original_train": 24,
            "original_validation": 16,
            "final_fit_count": 40,
            "final_test_count": 10,
        }
    )
    config["task"]["feature_count_before_filter"] = 6
    config["pipeline"]["gene_filter"]["required_fit_sample_count"] = 8
    return config


def synthetic_source(*, fail_stage: str | None = None) -> SyntheticFixtureSource:
    train_y = np.repeat(LABELS, 8)
    test_y = np.repeat(LABELS, 2)
    train = np.zeros((40, 6), dtype=float)
    train[:8, 0] = 2.0  # Exactly 8 > 1: retained at the 20% boundary.
    train[:7, 1] = 2.0  # Only 7 > 1: rejected.
    train[:, 2] = 1.0  # Equality does not pass the strict > 1 rule.
    train[:, 3] = 2.0 + np.repeat(np.arange(5), 8)
    train[:, 4] = np.random.default_rng(12).normal(2, 0.1, 40)
    train[:, 5] = 0.5
    test = np.zeros((10, 6), dtype=float)
    test[:, 3] = 2.0 + np.repeat(np.arange(5), 2)
    test[:, 4] = 2.0
    genes = [f"G{i}" for i in range(6)]
    fit_samples = tuple(f"SYN-FS-{i}" for i in range(40))
    test_samples = tuple(f"SYN-TS-{i}" for i in range(10))
    fit = FitData(
        pd.DataFrame(train, index=fit_samples, columns=genes),
        train_y,
        tuple(f"SYN-FP-{i}" for i in range(40)),
        fit_samples,
    )
    held = HeldOutData(
        pd.DataFrame(test, index=test_samples, columns=genes),
        tuple(f"SYN-TP-{i}" for i in range(10)),
        test_samples,
    )
    return SyntheticFixtureSource(fit, held, test_y, fail_stage=fail_stage)


def run_fixture(
    tmp_path: Path, source: SyntheticFixtureSource | None = None
) -> tuple[dict, SyntheticFixtureSource]:
    source = source or synthetic_source()
    manifest = execute_once(
        source,
        synthetic_config(),
        tmp_path / "result",
        config_sha256="synthetic-config-hash",
        source_revision="synthetic-revision",
        authorization_sha256=None,
    )
    return manifest, source


def states(path: Path) -> list[str]:
    return [json.loads(line)["state"] for line in path.read_text(encoding="utf-8").splitlines()]


def test_memberships_are_exact_and_disjoint() -> None:
    parts = {
        "train": (("SYN-P1", "SYN-S1"), ("SYN-P2", "SYN-S2")),
        "validation": (("SYN-P3", "SYN-S3"),),
        "test": (("SYN-P4", "SYN-S4"),),
    }
    fit, test = validate_memberships(parts, {"train": 2, "validation": 1, "test": 1})
    assert fit == parts["train"] + parts["validation"]
    assert test == parts["test"]
    bad = copy.deepcopy(parts)
    bad["test"] = (("SYN-P1", "SYN-S4"),)
    with pytest.raises(ValueError, match="overlapping"):
        validate_memberships(bad, {"train": 2, "validation": 1, "test": 1})
    with pytest.raises(ValueError, match="count"):
        validate_memberships(parts, {"train": 3, "validation": 1, "test": 1})


def test_config_and_threshold_boundaries() -> None:
    config = synthetic_config()
    validate_protocol(config, synthetic_fixture=True)
    with pytest.raises(ValueError, match="Real final-test population"):
        validate_protocol(config)
    config["pipeline"]["classifier"]["C"] = 0.1
    with pytest.raises(ValueError, match="approved scientific protocol"):
        validate_protocol(config, synthetic_fixture=True)


def test_one_shot_state_audit_metrics_probabilities_and_reproducibility(tmp_path: Path) -> None:
    first, source = run_fixture(tmp_path / "first")
    output = tmp_path / "first" / "result"
    assert states(output / "run_state.jsonl") == [
        "reserved",
        "fitting",
        "fitted",
        "test_access_started",
        "completed",
    ]
    assert source.events == [
        "preflight",
        "fit_loaded",
        "test_expression_loaded",
        "predictions_saved",
        "test_labels_loaded",
    ]
    fit = json.loads((output / "training_manifest.json").read_text(encoding="utf-8"))
    assert fit["required_prevalence_count"] == 8
    assert fit["fit_count"] == fit["scaler_fit_count"] == 40
    genes = [
        row["gene_id"]
        for row in csv.DictReader(
            (output / "training_gene_list.tsv").open(encoding="utf-8"), delimiter="\t"
        )
    ]
    assert "G0" in genes and "G1" not in genes and "G2" not in genes
    rows = list(
        csv.DictReader((output / "test_predictions.tsv").open(encoding="utf-8"), delimiter="\t")
    )
    assert len(rows) == len({row["patient_id"] for row in rows}) == 10
    assert list(rows[0]) == [
        "patient_id",
        "sample_id",
        "predicted_label",
        *(f"probability_{label}" for label in LABELS),
    ]
    for row in rows:
        assert sum(float(row[f"probability_{label}"]) for label in LABELS) == pytest.approx(1)
    metrics = json.loads((output / "test_metrics.json").read_text(encoding="utf-8"))
    matrix = np.asarray(metrics["confusion_matrix"])
    assert matrix.shape == (5, 5)
    assert matrix.sum() == 10
    assert matrix.sum(axis=1).tolist() == [2] * 5
    assert matrix.sum(axis=0).tolist() == [
        metrics["predicted_class_counts"][label] for label in LABELS
    ]
    assert first["source_revision"] == "synthetic-revision"
    assert first["config_sha256"] == "synthetic-config-hash"
    assert first["model_sha256"] == first["output_hashes"]["pipeline.joblib"]
    assert all((output / name).exists() for name in first["output_hashes"])
    with pytest.raises(FileExistsError):
        run_fixture(tmp_path / "first", source)
    second, _ = run_fixture(tmp_path / "second")
    output2 = tmp_path / "second" / "result"
    assert (output / "test_predictions.tsv").read_bytes() == (
        output2 / "test_predictions.tsv"
    ).read_bytes()
    assert (output / "test_metrics.json").read_bytes() == (
        output2 / "test_metrics.json"
    ).read_bytes()
    assert first["fit_identity_sha256"] == second["fit_identity_sha256"]
    assert first["ordered_gene_list_sha256"] == second["ordered_gene_list_sha256"]


def test_failure_before_test_access_is_preserved(tmp_path: Path) -> None:
    source = synthetic_source(fail_stage="fit_load")
    with pytest.raises(RuntimeError, match="fit-load"):
        run_fixture(tmp_path, source)
    output = tmp_path / "result"
    assert states(output / "run_state.jsonl") == ["reserved", "fitting", "failed"]
    assert source.events == ["preflight", "fit_loaded"]
    assert json.loads((output / "failure.json").read_text())["stage"] == "fitting"
    assert not (output / "test_predictions.tsv").exists()
    with pytest.raises(FileExistsError):
        run_fixture(tmp_path)


def test_preflight_failure_is_preserved_without_fit_or_test_access(tmp_path: Path) -> None:
    source = synthetic_source(fail_stage="preflight")
    with pytest.raises(RuntimeError, match="preflight"):
        run_fixture(tmp_path, source)
    output = tmp_path / "result"
    assert states(output / "run_state.jsonl") == ["reserved", "failed"]
    assert source.events == ["preflight"]
    assert json.loads((output / "failure.json").read_text())["stage"] == "preflight"
    with pytest.raises(FileExistsError):
        run_fixture(tmp_path)


def test_failure_after_test_access_is_preserved(tmp_path: Path) -> None:
    source = synthetic_source(fail_stage="test_labels")
    with pytest.raises(RuntimeError, match="test-label"):
        run_fixture(tmp_path, source)
    output = tmp_path / "result"
    assert states(output / "run_state.jsonl") == [
        "reserved",
        "fitting",
        "fitted",
        "test_access_started",
        "failed",
    ]
    assert (output / "test_predictions.tsv").exists()
    assert source.events[-2:] == ["predictions_saved", "test_labels_loaded"]
    assert json.loads((output / "failure.json").read_text())["stage"] == "test_access_started"
    with pytest.raises(FileExistsError):
        run_fixture(tmp_path)


def test_fixed_metric_orientation_and_zero_division() -> None:
    truth = np.asarray(LABELS, dtype=object)
    pred = np.asarray(["Luminal A"] * 5, dtype=object)
    result = _score(truth, pred, LABELS)
    assert result["macro_f1"] == pytest.approx(1 / 15)
    assert result["balanced_accuracy"] == pytest.approx(0.2)
    assert result["accuracy"] == pytest.approx(0.2)
    assert result["per_class"]["Normal-like"]["precision"] == 0
    assert result["confusion_matrix"] == [[1, 0, 0, 0, 0]] * 5
    assert result["row_normalized_confusion_matrix"] == [[1, 0, 0, 0, 0]] * 5


def test_real_execution_rejected_before_source_access(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    config_path = Path(__file__).resolve().parents[1] / "configs/final_test_v1.yaml"
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    with pytest.raises(AuthorizationError, match="No separate"):
        run_real(tmp_path, config_path)
    assert not (tmp_path / "data/processed/final_test_v1").exists()
    monkeypatch.setattr(RealSource, "preflight", lambda *a: pytest.fail("Real data accessed"))
    with pytest.raises(AuthorizationError, match="disabled"):
        execute_once(
            RealSource(tmp_path),
            config,
            tmp_path / "result",
            config_sha256="x",
            source_revision="x",
            authorization_sha256="fake",
        )
    assert not (tmp_path / "result").exists()
    with pytest.raises(ValueError, match="SYN-"):
        SyntheticFixtureSource(
            FitData(
                pd.DataFrame([[2.0]], index=["TCGA-real"], columns=["G"]),
                np.asarray([LABELS[0]]),
                ("TCGA-real",),
                ("TCGA-real",),
            ),
            HeldOutData(
                pd.DataFrame([[2.0]], index=["SYN-T"], columns=["G"]), ("SYN-T",), ("SYN-T",)
            ),
            np.asarray([LABELS[0]]),
        )
