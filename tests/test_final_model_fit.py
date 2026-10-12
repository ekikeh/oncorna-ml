"""Synthetic-only tests for the sealed fitting stage and locked real boundary."""

from __future__ import annotations

import copy
import hashlib
import json
import uuid
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest
import yaml

import oncorna.final_model_fit as fit_module
from oncorna.final_fit_export import sequence_digest
from oncorna.final_model_fit import (
    FitInputs,
    expected_authorization,
    run_real_fit,
    run_synthetic_fit,
    verify_authorization_record,
    verify_sealed_fit,
)
from oncorna.final_test import LABELS, AuthorizationError, run_real
from oncorna.preprocessing import identity_digest, sha256_file


def fixture(tmp_path: Path) -> tuple[FitInputs, dict, dict]:
    config = copy.deepcopy(
        yaml.safe_load(
            (Path(__file__).resolve().parents[1] / "configs/final_test_v1.yaml").read_bytes()
        )
    )
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
    values = np.zeros((40, 6), dtype="<f8")
    values[:8, 0] = 2.0  # Boundary: exactly 20% retained.
    values[:7, 1] = 2.0
    values[:, 2] = 1.0  # Equality does not pass strict > 1.
    values[:, 3] = 2.0 + np.repeat(np.arange(5), 8)
    values[:, 4] = np.random.default_rng(12).normal(2, 0.1, 40)
    values[:, 5] = 0.5
    path = tmp_path / "synthetic_fit.npy"
    np.save(path, values)
    pairs = tuple((f"SYN-P-{i}", f"SYN-S-{i}") for i in range(40))
    genes = tuple(f"SYN-G-{i}" for i in range(6))
    labels = tuple(np.repeat(LABELS, 8))
    inputs = FitInputs(
        path,
        pairs,
        genes,
        labels,
        "true_label",
        {"training_oof_predictions": "a" * 64, "validation_predictions": "b" * 64},
    )
    kwargs = {
        "artifact_sha256": sha256_file(path),
        "config_sha256": "c" * 64,
        "split_sha256": "d" * 64,
        "expected_fit_identity_sha256": identity_digest(pairs),
        "expected_gene_sha256": sequence_digest(genes),
        "expected_true_label_sha256": hashlib.sha256(
            json.dumps(labels, ensure_ascii=False, separators=(",", ":")).encode()
        ).hexdigest(),
    }
    return inputs, config, kwargs


def test_frozen_pipeline_and_complete_seal(tmp_path: Path) -> None:
    inputs, config, kwargs = fixture(tmp_path)
    output = tmp_path / "sealed_fit"
    manifest = run_synthetic_fit(inputs, config, output, **kwargs)
    assert verify_sealed_fit(output, manifest) == manifest
    assert manifest["fit_count"] == 40
    assert manifest["retained_gene_count"] == 3
    training = json.loads((output / "training_manifest.json").read_bytes())
    assert training["retained_gene_ids"] == ["SYN-G-0", "SYN-G-3", "SYN-G-4"]
    assert training["prevalence_required_count"] == 8
    assert training["prevalence_expression_gt"] == 1.0
    assert training["scaler_sample_count"] == 40
    assert training["fitted_sample_count"] == 40
    assert training["classifier_classes"] == sorted(LABELS)
    assert manifest["model_parameters"]["C"] == 1.0
    assert manifest["model_parameters"]["solver"] == "lbfgs"
    assert manifest["full_source_hash_verified"] is False
    assert manifest["source_row_equality_verified"] is False
    assert not (output / "test_predictions.tsv").exists()
    assert not (output / "test_metrics.json").exists()


@pytest.mark.parametrize("change", ["label_column", "patient", "gene", "label_digest"])
def test_wrong_labels_membership_or_gene_order_rejected(tmp_path: Path, change: str) -> None:
    inputs, config, kwargs = fixture(tmp_path)
    if change == "label_column":
        inputs = replace(inputs, label_column="predicted_label")
    elif change == "patient":
        inputs = replace(inputs, fit_pairs=tuple(reversed(inputs.fit_pairs)))
    elif change == "gene":
        inputs = replace(inputs, gene_ids=tuple(reversed(inputs.gene_ids)))
    else:
        kwargs["expected_true_label_sha256"] = "0" * 64
    with pytest.raises(ValueError):
        run_synthetic_fit(inputs, config, tmp_path / "rejected", **kwargs)
    assert not (tmp_path / "rejected").exists()


def test_changed_expression_and_wrong_artifact_digest_rejected(tmp_path: Path) -> None:
    inputs, config, kwargs = fixture(tmp_path)
    with inputs.expression_path.open("r+b") as handle:
        handle.seek(-8, 2)
        handle.write(b"changed!")
    with pytest.raises(ValueError, match="checksum"):
        run_synthetic_fit(inputs, config, tmp_path / "rejected", **kwargs)
    assert not (tmp_path / "rejected").exists()


def test_tampering_and_incomplete_state_are_not_success(tmp_path: Path) -> None:
    inputs, config, kwargs = fixture(tmp_path)
    output = tmp_path / "sealed_fit"
    run_synthetic_fit(inputs, config, output, **kwargs)
    with (output / "pipeline.joblib").open("ab") as handle:
        handle.write(b"tamper")
    with pytest.raises(ValueError, match="checksum"):
        verify_sealed_fit(output)
    with pytest.raises(FileExistsError):
        run_synthetic_fit(inputs, config, output, **kwargs)


def test_interrupted_attempt_preserved_without_retry(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    inputs, config, kwargs = fixture(tmp_path)
    output = tmp_path / "interrupted"
    original = fit_module._new_json

    def interrupt(path: Path, value: dict) -> None:
        if path.name == "success.json":
            raise OSError("synthetic interruption")
        original(path, value)

    monkeypatch.setattr(fit_module, "_new_json", interrupt)
    with pytest.raises(OSError, match="synthetic interruption"):
        run_synthetic_fit(inputs, config, output, **kwargs)
    assert (output / "failure.json").exists()
    with pytest.raises(ValueError, match="incomplete"):
        verify_sealed_fit(output)
    with pytest.raises(FileExistsError):
        run_synthetic_fit(inputs, config, output, **kwargs)


def test_success_marker_and_terminal_state_must_agree(tmp_path: Path) -> None:
    inputs, config, kwargs = fixture(tmp_path)
    output = tmp_path / "sealed_fit"
    run_synthetic_fit(inputs, config, output, **kwargs)
    (output / "success.json").write_text('{"status":"sealed","manifest_sha256":"bad"}')
    with pytest.raises(ValueError, match="success marker"):
        verify_sealed_fit(output)


def test_training_manifest_and_audit_history_tampering_rejected(tmp_path: Path) -> None:
    inputs, config, kwargs = fixture(tmp_path)
    output = tmp_path / "sealed_fit"
    run_synthetic_fit(inputs, config, output, **kwargs)
    training_path = output / "training_manifest.json"
    training = json.loads(training_path.read_bytes())
    training["retained_gene_ids"] = ["SYN-G-1"]
    training_path.write_text(json.dumps(training))
    with pytest.raises(ValueError, match="checksum"):
        verify_sealed_fit(output)
    training_path.write_text("{}")
    (output / "fit_state.jsonl").write_text('{"state":"completed"}\n')
    with pytest.raises(ValueError, match="incomplete"):
        verify_sealed_fit(output)


def test_authorization_contract_rejects_revision_type_path_and_dirty_checkout(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    (tmp_path / "data/processed").mkdir(parents=True)
    record = expected_authorization(
        tmp_path,
        revision="a" * 40,
        attempt_id=str(uuid.uuid4()),
        approval_reference="synthetic human approval",
        config_sha256="b" * 64,
        split_sha256="c" * 64,
        fit_expression_sha256="d" * 64,
        ordered_fit_identity_sha256="e" * 64,
        ordered_gene_list_sha256="f" * 64,
        ordered_true_label_sha256="1" * 64,
        label_table_hashes={
            "training_oof_predictions": "2" * 64,
            "validation_predictions": "3" * 64,
        },
    )
    verify_authorization_record(record, record)
    for altered in (
        {**record, "approved_commit": "0" * 40},
        {**record, "canonical_output": str(tmp_path / "wrong")},
        {**record, "model_parameters": {**record["model_parameters"], "fit_intercept": 1}},
        {**record, "extra": "unexpected"},
    ):
        with pytest.raises(AuthorizationError):
            verify_authorization_record(altered, record)

    def fake_git(args: list[str], **_kwargs: object) -> str:
        if args[-1] == "--show-toplevel":
            return str(tmp_path) + "\n"
        if args[-1] == "HEAD":
            return "a" * 40 + "\n"
        return "dirty\n"

    monkeypatch.setattr(fit_module.subprocess, "check_output", fake_git)
    auth_path = tmp_path / fit_module.AUTHORIZATION_PATH
    auth_path.write_text(json.dumps(record))
    with pytest.raises(AuthorizationError):
        fit_module.verify_real_fit_authorization(tmp_path, auth_path, record)


def test_prior_output_or_claim_blocks_authorization(tmp_path: Path) -> None:
    (tmp_path / "data/processed").mkdir(parents=True)
    auth = tmp_path / fit_module.AUTHORIZATION_PATH
    auth.write_text("{}")
    output = tmp_path / fit_module.OUTPUT_DIRECTORY
    output.mkdir()
    with pytest.raises(FileExistsError):
        fit_module.verify_real_fit_authorization(tmp_path, auth, {})
    output.rmdir()
    (tmp_path / fit_module.ATTEMPT_PATH).write_text("{}")
    with pytest.raises(FileExistsError):
        fit_module.verify_real_fit_authorization(tmp_path, auth, {})


def test_real_fit_and_final_test_remain_locked(tmp_path: Path) -> None:
    with pytest.raises(AuthorizationError, match="disabled"):
        run_real_fit(tmp_path)
    config_path = Path(__file__).resolve().parents[1] / "configs/final_test_v1.yaml"
    with pytest.raises(AuthorizationError):
        run_real(tmp_path, config_path)
    assert not (tmp_path / fit_module.OUTPUT_DIRECTORY).exists()
