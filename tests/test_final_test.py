"""Synthetic-only checks for the frozen final-test infrastructure."""

from __future__ import annotations

import copy
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import yaml

import oncorna.final_test as final_test_module
from oncorna.final_fit_export import ExportPlan, export_synthetic
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


def synthetic_real_adapter_metadata(root: Path) -> tuple[dict, SyntheticFixtureSource]:
    """Build metadata and fit-only labels; deliberately create no full-cohort value files."""
    config = synthetic_config()
    fixture = synthetic_source()
    train_pairs = list(zip(fixture.fit.patient_ids[:24], fixture.fit.sample_ids[:24], strict=True))
    validation_pairs = list(
        zip(fixture.fit.patient_ids[24:], fixture.fit.sample_ids[24:], strict=True)
    )
    test_pairs = list(zip(fixture.test.patient_ids, fixture.test.sample_ids, strict=True))

    def write(path: Path, content: str) -> str:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content, encoding="utf-8", newline="\n")
        return hashlib.sha256(path.read_bytes()).hexdigest()

    default_hash = write(root / "configs/default.yaml", "synthetic: true\n")
    matrix_manifest_hash = write(root / "data/processed/ml_matrix_manifest.json", "{}\n")
    split = {
        "split_id": "split_v1",
        "split_unit": "patient_id",
        "configuration": {"sha256": default_hash},
        "expression_matrix_header_validation": {"matrix_manifest_sha256": matrix_manifest_hash},
        "partitions": {
            name: {
                "n_patients": len(pairs),
                "patient_ids": [patient for patient, _ in pairs],
                "sample_ids": [sample for _, sample in pairs],
                "patient_sample_map": [
                    {"patient_id": patient, "sample_id": sample} for patient, sample in pairs
                ],
            }
            for name, pairs in (
                ("train", train_pairs),
                ("validation", validation_pairs),
                ("test", test_pairs),
            )
        },
    }
    schema = {
        "orientation": "samples_by_genes",
        "expression_units": config["task"]["expression_units"],
        "shape": [50, 6],
        "dtype": "float64",
        "gene_ids": list(fixture.fit.X.columns),
        "sample_ids": list(fixture.fit.sample_ids + fixture.test.sample_ids),
    }
    oof = root / "data/processed/modeling_v1/oof_predictions.tsv"
    validation = root / "data/processed/validation_v1/validation_predictions.tsv"
    write(
        oof,
        "patient_id\tsample_id\ttrue_label\tmodel\tC\n"
        + "".join(
            f"{patient}\t{sample}\t{label}\tlogistic_l2_multinomial\t1.0\n"
            for (patient, sample), label in zip(train_pairs, fixture.fit.y[:24], strict=True)
        ),
    )
    write(
        validation,
        "patient_id\tsample_id\ttrue_label\n"
        + "".join(
            f"{patient}\t{sample}\t{label}\n"
            for (patient, sample), label in zip(validation_pairs, fixture.fit.y[24:], strict=True)
        ),
    )
    manifests = {
        "phase2b_manifest": {
            "output_hashes": {oof.name: hashlib.sha256(oof.read_bytes()).hexdigest()}
        },
        "phase2c_manifest": {
            "output_hashes": {validation.name: hashlib.sha256(validation.read_bytes()).hexdigest()}
        },
    }
    for key, path in config["inputs"].items():
        if key == "expected_sha256" or key in final_test_module.TEST_CONTAINING_INPUTS:
            continue
        target = root / path
        if key == "frozen_split":
            content = json.dumps(split)
        elif key == "all_gene_schema":
            content = json.dumps(schema)
        elif key in manifests:
            content = json.dumps(manifests[key])
        else:
            content = "synthetic metadata\n"
        config["inputs"]["expected_sha256"][key] = write(target, content)
    synthetic_array = root / "synthetic_source.npy"
    np.save(synthetic_array, np.vstack([fixture.fit.X.to_numpy(), fixture.test.X.to_numpy()]))
    export_synthetic(
        ExportPlan(
            source_array=synthetic_array,
            output_dir=root / final_test_module.FIT_EXPRESSION_PATH.parent,
            fit_pairs=tuple(train_pairs + validation_pairs),
            test_pairs=tuple(test_pairs),
            sample_ids=tuple(schema["sample_ids"]),
            gene_ids=tuple(schema["gene_ids"]),
            frozen_split_sha256=config["inputs"]["expected_sha256"]["frozen_split"],
            source_schema_sha256=config["inputs"]["expected_sha256"]["all_gene_schema"],
            source_all_gene_expression_sha256=config["inputs"]["expected_sha256"][
                "all_gene_expression"
            ],
            source_expression_sha256=config["inputs"]["expected_sha256"]["source_expression"],
            approved_cohort_identifier="TCGA-BRCA Xena primary-tumor cohort",
            approved_cohort_path=config["inputs"]["cohort"],
            approved_cohort_sha256=config["inputs"]["expected_sha256"]["cohort"],
            expression_units=config["task"]["expression_units"],
            command="synthetic fixture export",
        )
    )
    return config, fixture


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
        "test_access_authorized",
        "full_provenance_verified",
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
    assert (output / "success.json").exists()
    assert first["output_hashes"].get("run_state.jsonl") is None
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
    completed = json.loads((output / "run_state.jsonl").read_text().splitlines()[-1])
    assert completed["manifest_sha256"] == final_test_module.sha256_file(output / "manifest.json")
    assert completed["success_sha256"] == final_test_module.sha256_file(output / "success.json")
    assert (
        json.loads((output / "success.json").read_text())["manifest_sha256"]
        == completed["manifest_sha256"]
    )
    for name, digest in first["output_hashes"].items():
        assert final_test_module.sha256_file(output / name) == digest
    assert final_test_module.verify_completed_run(output) == first


def test_held_out_values_cannot_change_fitted_artifacts(tmp_path: Path) -> None:
    baseline, _ = run_fixture(tmp_path / "baseline")
    changed_source = synthetic_source()
    changed_source.test.X.iloc[:, :] = 9999.0
    changed_source.test_labels = np.roll(changed_source.test_labels, 2)
    changed, _ = run_fixture(tmp_path / "changed", changed_source)
    assert baseline["fit_identity_sha256"] == changed["fit_identity_sha256"]
    assert baseline["ordered_gene_list_sha256"] == changed["ordered_gene_list_sha256"]
    assert baseline["model_sha256"] == changed["model_sha256"]


def test_real_adapter_preflight_and_fit_labels_use_only_safe_artifacts(tmp_path: Path) -> None:
    config, fixture = synthetic_real_adapter_metadata(tmp_path)
    source = RealSource(tmp_path)
    assert not (tmp_path / config["inputs"]["cohort"]).exists()
    assert not (tmp_path / config["inputs"]["all_gene_expression"]).exists()
    snapshot = source.preflight(config)
    assert len(snapshot.fit_pairs) == 40
    assert len(snapshot.test_pairs) == 10
    assert len(snapshot.gene_ids) == 6
    assert not set(snapshot.input_hashes) & final_test_module.TEST_CONTAINING_INPUTS
    fit = source.load_fit(snapshot)
    assert fit.y.tolist() == fixture.fit.y.tolist()
    assert np.array_equal(fit.X.to_numpy(), fixture.fit.X.to_numpy())
    with pytest.raises(AuthorizationError, match="test-access stage"):
        source.verify_full_provenance(snapshot)
    with pytest.raises(AuthorizationError, match="test-access stage"):
        source.load_test_expression(snapshot)
    with pytest.raises(AuthorizationError, match="sealed predictions"):
        source.load_test_labels(snapshot)
    source.begin_test_access()
    with pytest.raises(FileNotFoundError):
        source.verify_full_provenance(snapshot)
    assert not source.full_provenance_verified
    with pytest.raises(AuthorizationError, match="test-access stage"):
        source.load_test_expression(snapshot)


def test_fit_only_artifact_change_is_rejected_before_fitting(tmp_path: Path) -> None:
    config, _ = synthetic_real_adapter_metadata(tmp_path)
    source = RealSource(tmp_path)
    snapshot = source.preflight(config)
    path = tmp_path / final_test_module.FIT_EXPRESSION_PATH
    with path.open("ab") as handle:
        handle.write(b"tampered")
    with pytest.raises(ValueError, match="changed after preflight"):
        source.load_fit(snapshot)
    assert not source.test_access_started


@pytest.mark.parametrize("interruption", ["before_manifest", "after_manifest"])
def test_interrupted_sealing_never_records_completion(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, interruption: str
) -> None:
    original = final_test_module._atomic_write_json_new

    def interrupt(path: Path, value: dict) -> None:
        if (interruption == "before_manifest" and path.name == "manifest.json") or (
            interruption == "after_manifest" and path.name == "success.json"
        ):
            raise KeyboardInterrupt("synthetic crash")
        original(path, value)

    monkeypatch.setattr(final_test_module, "_atomic_write_json_new", interrupt)
    with pytest.raises(KeyboardInterrupt, match="synthetic crash"):
        run_fixture(tmp_path)
    output = tmp_path / "result"
    assert states(output / "run_state.jsonl")[-1] == "test_access_started"
    assert (output / "manifest.json").exists() == (interruption == "after_manifest")
    assert not (output / "success.json").exists()
    with pytest.raises(ValueError, match="no terminal completed"):
        final_test_module.verify_completed_run(output)
    with pytest.raises(FileExistsError):
        run_fixture(tmp_path)


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
