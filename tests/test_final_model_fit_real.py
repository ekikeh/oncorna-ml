"""Synthetic-only coverage of real fitting metadata and authorization boundaries."""

from __future__ import annotations

import json
import os
import uuid
from pathlib import Path

import numpy as np
import pytest
import yaml
from test_final_test import synthetic_real_adapter_metadata

import oncorna.final_model_fit as fit_module
import oncorna.final_model_fit_real as real
from oncorna.final_fit_export import sequence_digest
from oncorna.final_model_fit import (
    ATTEMPT_PATH,
    AUTHORIZATION_PATH,
    OUTPUT_DIRECTORY,
    _digest,
    expected_authorization,
    run_real_fit,
    verify_sealed_fit,
)
from oncorna.final_test import LABELS, AuthorizationError, run_real
from oncorna.preprocessing import identity_digest, sha256_file


def synthetic_binding(tmp_path: Path) -> tuple[dict, real.VerifiedFitBinding]:
    config, _ = synthetic_real_adapter_metadata(tmp_path)
    path = tmp_path / "configs/final_test_v1.yaml"
    path.write_text(yaml.safe_dump(config), encoding="utf-8")
    counts = dict.fromkeys(LABELS, 8)
    binding = real._derive_binding(tmp_path, config, synthetic_fixture=True, class_counts=counts)
    return config, binding


def synthetic_approval(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> tuple[Path, real.VerifiedFitBinding]:
    _, binding = synthetic_binding(tmp_path)
    monkeypatch.setattr(real, "derive_real_fit_binding", lambda _root: binding)
    monkeypatch.setattr(real, "_revision", lambda _root: "a" * 40)
    authorization = tmp_path / AUTHORIZATION_PATH
    record = expected_authorization(
        tmp_path,
        revision="a" * 40,
        attempt_id=str(uuid.uuid4()),
        approval_reference="synthetic explicit approval",
        config_sha256=binding.config_sha256,
        split_sha256=binding.split_sha256,
        fit_expression_sha256=binding.expression_sha256,
        ordered_fit_identity_sha256=identity_digest(binding.inputs.fit_pairs),
        ordered_gene_list_sha256=sequence_digest(binding.inputs.gene_ids),
        ordered_true_label_sha256=_digest(binding.inputs.true_labels),
        label_table_hashes=binding.inputs.label_table_hashes,
    )
    authorization.write_text(json.dumps(record), encoding="utf-8")
    return authorization, binding


def test_frozen_metadata_and_true_label_derivation(tmp_path: Path) -> None:
    config, binding = synthetic_binding(tmp_path)
    assert len(binding.inputs.fit_pairs) == 40
    assert binding.test_count == 10
    assert len(binding.inputs.gene_ids) == 6
    assert binding.inputs.true_labels == tuple(label for label in LABELS for _ in range(8))
    assert binding.inputs.label_column == "true_label"
    assert binding.split_sha256 == config["inputs"]["expected_sha256"]["frozen_split"]
    assert binding.expression_sha256 == sha256_file(binding.inputs.expression_path)
    assert not (tmp_path / config["inputs"]["cohort"]).exists()
    assert not (tmp_path / config["inputs"]["all_gene_expression"]).exists()


def test_predictions_are_not_fitting_targets(tmp_path: Path) -> None:
    config, binding = synthetic_binding(tmp_path)
    path = tmp_path / "data/processed/modeling_v1/oof_predictions.tsv"
    text = path.read_text(encoding="utf-8")
    path.write_text(
        text.replace("\ttrue_label\t", "\ttrue_label\tpredicted_label\t").replace(
            "\tlogistic_l2_multinomial", "\twrong_prediction\tlogistic_l2_multinomial"
        ),
        encoding="utf-8",
    )
    # A changed table is rejected until its training-side provenance is updated.
    with pytest.raises(ValueError, match="checksum"):
        real._derive_binding(
            tmp_path, config, synthetic_fixture=True, class_counts=dict.fromkeys(LABELS, 8)
        )
    manifest_path = tmp_path / config["inputs"]["phase2b_manifest"]
    manifest = json.loads(manifest_path.read_text())
    manifest["output_hashes"][path.name] = sha256_file(path)
    manifest_path.write_text(json.dumps(manifest))
    config["inputs"]["expected_sha256"]["phase2b_manifest"] = sha256_file(manifest_path)
    updated = real._derive_binding(
        tmp_path, config, synthetic_fixture=True, class_counts=dict.fromkeys(LABELS, 8)
    )
    assert updated.inputs.true_labels == binding.inputs.true_labels
    assert updated.inputs.true_labels[0] == LABELS[0]


@pytest.mark.parametrize("change", ["membership", "gene", "counts", "expression"])
def test_modified_frozen_inputs_rejected(tmp_path: Path, change: str) -> None:
    config, binding = synthetic_binding(tmp_path)
    if change in {"membership", "gene"}:
        key = "frozen_split" if change == "membership" else "all_gene_schema"
        path = tmp_path / config["inputs"][key]
        data = json.loads(path.read_text())
        if change == "membership":
            data["partitions"]["train"]["patient_sample_map"][0]["patient_id"] = "WRONG"
        else:
            data["gene_ids"] = list(reversed(data["gene_ids"]))
        path.write_text(json.dumps(data))
    elif change == "counts":
        config["partitions"]["original_train"] = 23
    else:
        with binding.inputs.expression_path.open("ab") as handle:
            handle.write(b"changed")
    with pytest.raises((ValueError, KeyError)):
        real._derive_binding(
            tmp_path, config, synthetic_fixture=True, class_counts=dict.fromkeys(LABELS, 8)
        )


def test_authorization_is_derived_and_claim_is_exclusive(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    authorization, binding = synthetic_approval(tmp_path, monkeypatch)
    prepared = real.verify_real_fit_preparation(tmp_path, authorization)
    assert prepared.binding == binding
    assert prepared.authorization_sha256 == sha256_file(authorization)
    with pytest.raises(AuthorizationError, match="redirected"):
        real.verify_real_fit_preparation(tmp_path, tmp_path / "wrong.json")
    real.claim_real_fit_attempt(tmp_path, authorization)
    claim = json.loads((tmp_path / ATTEMPT_PATH).read_text())
    assert claim["authorization_sha256"] == sha256_file(authorization)
    assert claim["status"] == "claimed_before_fitting"
    assert not (tmp_path / OUTPUT_DIRECTORY).exists()
    with pytest.raises(FileExistsError):
        real.claim_real_fit_attempt(tmp_path, authorization)


@pytest.mark.parametrize(
    "field", ["approved_commit", "ordered_true_label_sha256", "canonical_output"]
)
def test_mismatched_authorization_rejected(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, field: str
) -> None:
    authorization, _ = synthetic_approval(tmp_path, monkeypatch)
    record = json.loads(authorization.read_text())
    record[field] = "wrong"
    authorization.write_text(json.dumps(record))
    with pytest.raises(AuthorizationError):
        real.verify_real_fit_preparation(tmp_path, authorization)
    assert not (tmp_path / ATTEMPT_PATH).exists()


def test_dirty_checkout_and_claim_failure_stop_before_output(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    authorization, _ = synthetic_approval(tmp_path, monkeypatch)
    monkeypatch.setattr(
        real, "_revision", lambda _root: (_ for _ in ()).throw(AuthorizationError("dirty checkout"))
    )
    with pytest.raises(AuthorizationError, match="dirty"):
        real.claim_real_fit_attempt(tmp_path, authorization)
    monkeypatch.setattr(real, "_revision", lambda _root: "a" * 40)
    monkeypatch.setattr(
        real,
        "_create_local_claim",
        lambda *_args: (_ for _ in ()).throw(PermissionError("synthetic ACL failure")),
    )
    with pytest.raises(PermissionError, match="ACL"):
        real.claim_real_fit_attempt(tmp_path, authorization)
    assert not (tmp_path / OUTPUT_DIRECTORY).exists()


def test_interrupted_claim_is_preserved_without_retry(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    authorization, _ = synthetic_approval(tmp_path, monkeypatch)
    original = real._create_local_claim

    def interrupt(marker: Path, record: dict[str, str]) -> None:
        original(marker, record)
        raise OSError("synthetic interruption after claim")

    monkeypatch.setattr(real, "_create_local_claim", interrupt)
    with pytest.raises(OSError, match="interruption"):
        real.claim_real_fit_attempt(tmp_path, authorization)
    assert (tmp_path / ATTEMPT_PATH).exists()
    with pytest.raises(FileExistsError):
        real.claim_real_fit_attempt(tmp_path, authorization)
    assert not (tmp_path / OUTPUT_DIRECTORY).exists()


def test_real_fit_and_final_test_still_locked(tmp_path: Path) -> None:
    with pytest.raises(AuthorizationError, match="disabled"):
        run_real_fit(tmp_path)
    config_path = Path(__file__).resolve().parents[1] / "configs/final_test_v1.yaml"
    with pytest.raises(AuthorizationError):
        run_real(tmp_path, config_path)


def test_claimed_synthetic_fit_reuses_sealed_pipeline(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    authorization, binding = synthetic_approval(tmp_path, monkeypatch)
    manifest = real._run_claimed_fit(tmp_path, authorization, synthetic_fixture=True)
    output = tmp_path / OUTPUT_DIRECTORY
    assert manifest == verify_sealed_fit(output)
    assert manifest["fit_count"] == len(binding.inputs.fit_pairs) == 40
    assert manifest["original_gene_count"] == 6
    assert manifest["attempt_id"] == json.loads(authorization.read_text())["attempt_id"]
    assert manifest["scope"] == "claimed_synthetic"
    assert manifest["authorization_sha256"] == sha256_file(authorization)
    assert (tmp_path / ATTEMPT_PATH).exists()
    assert not (output / "test_predictions.tsv").exists()
    with pytest.raises(FileExistsError):
        real._run_claimed_fit(tmp_path, authorization, synthetic_fixture=True)


@pytest.mark.parametrize("change", ["missing_claim", "attempt_id", "auth_hash", "values"])
def test_direct_claimed_engine_rechecks_evidence_without_callbacks(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, change: str
) -> None:
    authorization, binding = synthetic_approval(tmp_path, monkeypatch)
    record = json.loads(authorization.read_text())
    if change != "missing_claim":
        real.claim_real_fit_attempt(tmp_path, authorization)
    values = np.load(binding.inputs.expression_path, allow_pickle=False)
    if change == "values":
        values = values.copy()
        values[0, 0] += 1.0
    with pytest.raises(AuthorizationError):
        fit_module._fit_and_seal(
            binding.inputs,
            yaml.safe_load((tmp_path / "configs/final_test_v1.yaml").read_bytes()),
            tmp_path / OUTPUT_DIRECTORY,
            values=values,
            artifact_sha256=binding.expression_sha256,
            config_sha256=binding.config_sha256,
            split_sha256=binding.split_sha256,
            source_revision="a" * 40,
            attempt_id=str(uuid.uuid4()) if change == "attempt_id" else record["attempt_id"],
            authorization_reference=record["approval_reference"],
            execution_scope="claimed_synthetic",
            authorization_sha256="f" * 64 if change == "auth_hash" else sha256_file(authorization),
            claimed_root=tmp_path,
            checkpoint=lambda _stage: None,
            validate_manifest=lambda _manifest: None,
        )
    assert not (tmp_path / OUTPUT_DIRECTORY).exists()


@pytest.mark.parametrize(
    "change",
    ["expression", "label_table", "claim", "replacement", "authorization", "split", "gene"],
)
def test_change_after_claim_consumes_attempt_before_fitting(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, change: str
) -> None:
    authorization, binding = synthetic_approval(tmp_path, monkeypatch)
    config = yaml.safe_load((tmp_path / "configs/final_test_v1.yaml").read_text())
    monkeypatch.setattr(
        real,
        "derive_real_fit_binding",
        lambda root: real._derive_binding(
            root, config, synthetic_fixture=True, class_counts=dict.fromkeys(LABELS, 8)
        ),
    )
    original_claim = real.claim_real_fit_attempt

    def change_after_claim(root: Path, path: Path) -> real.VerifiedPreparation:
        prepared = original_claim(root, path)
        if change == "expression":
            with binding.inputs.expression_path.open("ab") as handle:
                handle.write(b"mutation")
        elif change == "label_table":
            label_path = root / "data/processed/modeling_v1/oof_predictions.tsv"
            with label_path.open("a", encoding="utf-8") as handle:
                handle.write("tamper\n")
        elif change == "claim":
            (root / ATTEMPT_PATH).write_text("{}")
        elif change == "replacement":
            source = binding.inputs.expression_path
            replacement = source.with_name("replacement.npy")
            replacement.write_bytes(source.read_bytes())
            os.replace(replacement, source)
        elif change == "authorization":
            path.write_text(path.read_text() + " ")
        else:
            name = "frozen_split" if change == "split" else "all_gene_schema"
            metadata = root / config["inputs"][name]
            metadata.write_text(metadata.read_text() + " ")
        return prepared

    monkeypatch.setattr(real, "claim_real_fit_attempt", change_after_claim)
    with pytest.raises((AuthorizationError, ValueError)):
        real._run_claimed_fit(tmp_path, authorization, synthetic_fixture=True)
    assert (tmp_path / ATTEMPT_PATH).exists()
    assert not (tmp_path / OUTPUT_DIRECTORY).exists()
    with pytest.raises(FileExistsError):
        real._run_claimed_fit(tmp_path, authorization, synthetic_fixture=True)


def test_mutation_during_fit_preserves_failed_attempt(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    authorization, _ = synthetic_approval(tmp_path, monkeypatch)
    config = yaml.safe_load((tmp_path / "configs/final_test_v1.yaml").read_text())
    monkeypatch.setattr(
        real,
        "derive_real_fit_binding",
        lambda root: real._derive_binding(
            root, config, synthetic_fixture=True, class_counts=dict.fromkeys(LABELS, 8)
        ),
    )
    original_pipeline = fit_module.make_pipeline

    def mutate_after_fit(*args: object, **kwargs: object) -> object:
        pipeline = original_pipeline(*args, **kwargs)
        original_fit = pipeline.fit

        def fit_then_change(*fit_args: object, **fit_kwargs: object) -> object:
            result = original_fit(*fit_args, **fit_kwargs)
            label_path = tmp_path / "data/processed/validation_v1/validation_predictions.tsv"
            with label_path.open("a", encoding="utf-8") as handle:
                handle.write("mutation\n")
            return result

        pipeline.fit = fit_then_change
        return pipeline

    monkeypatch.setattr(fit_module, "make_pipeline", mutate_after_fit)
    with pytest.raises(ValueError, match="checksum"):
        real._run_claimed_fit(tmp_path, authorization, synthetic_fixture=True)
    output = tmp_path / OUTPUT_DIRECTORY
    assert (tmp_path / ATTEMPT_PATH).exists()
    assert (output / "failure.json").exists()
    with pytest.raises(ValueError, match="incomplete"):
        verify_sealed_fit(output)


def test_interrupted_seal_retains_claim_and_failure(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    authorization, _ = synthetic_approval(tmp_path, monkeypatch)
    original_write = fit_module._new_json

    def stop_at_success(path: Path, value: dict) -> None:
        if path.name == "success.json":
            raise OSError("synthetic seal interruption")
        original_write(path, value)

    monkeypatch.setattr(fit_module, "_new_json", stop_at_success)
    with pytest.raises(OSError, match="seal interruption"):
        real._run_claimed_fit(tmp_path, authorization, synthetic_fixture=True)
    assert (tmp_path / ATTEMPT_PATH).exists()
    assert (tmp_path / OUTPUT_DIRECTORY / "failure.json").exists()
    with pytest.raises(FileExistsError):
        real._run_claimed_fit(tmp_path, authorization, synthetic_fixture=True)


@pytest.mark.parametrize("name", ["pipeline.joblib", "manifest.json"])
def test_claimed_fit_rejects_tampered_seal(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, name: str
) -> None:
    authorization, _ = synthetic_approval(tmp_path, monkeypatch)
    real._run_claimed_fit(tmp_path, authorization, synthetic_fixture=True)
    output = tmp_path / OUTPUT_DIRECTORY
    with (output / name).open("ab") as handle:
        handle.write(b"tamper")
    with pytest.raises(ValueError):
        verify_sealed_fit(output)


@pytest.mark.parametrize("field", ["authorization_sha256", "attempt_id", "scope"])
def test_resealed_forged_claim_metadata_is_not_trusted(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, field: str
) -> None:
    authorization, _ = synthetic_approval(tmp_path, monkeypatch)
    real._run_claimed_fit(tmp_path, authorization, synthetic_fixture=True)
    output = tmp_path / OUTPUT_DIRECTORY
    manifest_path = output / "manifest.json"
    manifest = json.loads(manifest_path.read_bytes())
    manifest[field] = "f" * 64 if field == "authorization_sha256" else "forged"
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    manifest_sha = sha256_file(manifest_path)
    success_path = output / "success.json"
    success_path.write_text(
        json.dumps({"status": "sealed", "manifest_sha256": manifest_sha}), encoding="utf-8"
    )
    history_path = output / "fit_state.jsonl"
    history = [json.loads(line) for line in history_path.read_text().splitlines()]
    history[-1]["manifest_sha256"] = manifest_sha
    history[-1]["success_sha256"] = sha256_file(success_path)
    history_path.write_text("\n".join(json.dumps(row) for row in history) + "\n", encoding="utf-8")
    with pytest.raises((AuthorizationError, ValueError)):
        verify_sealed_fit(output)
