"""Synthetic checks for test-only model consumption and irreversible ordering."""

from __future__ import annotations

import copy
import json
import uuid
from pathlib import Path

import joblib
import numpy as np
import pytest
import yaml
from test_final_test import synthetic_config, synthetic_source

import oncorna.final_test_consumer as consumer
from oncorna.final_fit_export import inspect_npy_header_handle
from oncorna.final_test import (
    AuthorizationError,
    RealSource,
    SyntheticFileSource,
    _authenticated_source_rows,
    modeling_parameters,
)
from oncorna.final_test import TestData as HeldOutData
from oncorna.final_test_consumer import (
    _claim_payload,
    execute_test_only,
    expected_test_authorization,
    verify_approved_fit,
    verify_completed_test_only,
    verify_test_authorization,
)
from oncorna.modeling import make_pipeline
from oncorna.preprocessing import sha256_file


def prepared(tmp_path: Path):
    source = synthetic_source()
    config = synthetic_config()
    pipeline = make_pipeline(modeling_parameters(config), 1.0)
    pipeline.fit(source.fit.X, source.fit.y)  # Synthetic preparation only.
    path = tmp_path / "synthetic.joblib"
    joblib.dump(pipeline, path)
    genes = list(pipeline.named_steps["gene_filter"].selected_gene_ids_)
    fit_manifest = {
        "retained_gene_ids": genes,
        "fit_count": len(source.fit.y),
        "prevalence_required_count": 8,
        "classifier_classes": list(pipeline.named_steps["classifier"].classes_),
    }
    return source, config, path, fit_manifest


def consume(tmp_path: Path, source=None, *, loader=joblib.load):
    fixture, config, path, fit_manifest = prepared(tmp_path)
    if source is None:
        source = fixture
    result = execute_test_only(
        source,
        config,
        tmp_path / "result",
        path,
        fit_manifest,
        authorization_sha256="synthetic-authorization",
        model_sha256=sha256_file(path),
        source_revision="synthetic-revision",
        synthetic_model_loader=loader,
    )
    return result, source


def test_prediction_seal_precedes_label_access_and_reuse_is_rejected(tmp_path: Path):
    result, source = consume(tmp_path)
    assert source.events.index("predictions_saved") < source.events.index("test_labels_loaded")
    assert source.events.index("label_provenance_verified") < source.events.index(
        "test_labels_loaded"
    )
    assert result["test_count"] == 10
    entries = [
        json.loads(x) for x in (tmp_path / "result/run_state.jsonl").read_text().splitlines()
    ]
    assert [e["state"] for e in entries] == [
        "reserved",
        "test_access_started",
        "predictions_sealed",
        "completed",
    ]
    assert entries[2]["prediction_sha256"] == result["sealed_prediction_sha256"]
    assert result["class_order"] == [
        "Luminal A",
        "Luminal B",
        "Basal-like",
        "HER2-enriched",
        "Normal-like",
    ]
    with pytest.raises(FileExistsError):
        consume(tmp_path)


def test_file_backed_synthetic_labels_open_only_after_seal(tmp_path: Path):
    fixture, config, path, fit_manifest = prepared(tmp_path)
    labels = tmp_path / "synthetic_labels.tsv"
    labels.write_text(
        "patient_id\tsample_id\ttrue_label\n"
        + "".join(
            f"{patient}\t{sample}\t{label}\n"
            for patient, sample, label in zip(
                fixture.test.patient_ids, fixture.test.sample_ids, fixture.test_labels, strict=True
            )
        ),
        encoding="utf-8",
    )
    source = SyntheticFileSource(fixture.fit, fixture.test, labels)
    snapshot = source.preflight(config)
    with pytest.raises(RuntimeError, match="before predictions"):
        source.load_test_labels(snapshot)
    execute_test_only(
        source,
        config,
        tmp_path / "result",
        path,
        fit_manifest,
        authorization_sha256="synthetic-authorization",
        model_sha256=sha256_file(path),
        source_revision="synthetic-revision",
    )
    assert source.events.index("predictions_sealed_verified") < source.events.index(
        "label_file_opened"
    )
    assert source.events.index("label_provenance_verified") < source.events.index(
        "label_file_opened"
    )


def test_authenticated_source_uses_one_snapshot_and_rejects_changes(tmp_path: Path):
    fixture = synthetic_source()
    snapshot = fixture.preflight(synthetic_config())
    source = tmp_path / "source.npy"
    fit = tmp_path / "fit.npy"
    values = np.vstack((fixture.fit.X.to_numpy(), fixture.test.X.to_numpy()))
    np.save(source, values)
    np.save(fit, fixture.fit.X.to_numpy())
    with source.open("rb") as handle:
        offset = inspect_npy_header_handle(handle, values.shape).data_offset
    selected, stamp = _authenticated_source_rows(source, fit, snapshot, sha256_file(source), offset)
    assert np.array_equal(selected, fixture.test.X.to_numpy())
    assert not selected.flags.writeable
    assert stamp[2] == source.stat().st_size
    with pytest.raises(ValueError, match="checksum"):
        _authenticated_source_rows(source, fit, snapshot, "0" * 64, offset)
    changed = values.copy()
    changed[0, 0] += 1
    np.save(source, changed)
    with pytest.raises(ValueError, match="Fitting-only expression differs"):
        _authenticated_source_rows(source, fit, snapshot, sha256_file(source), offset)


@pytest.mark.parametrize("change", ["in_place", "replacement"])
def test_authenticated_test_rows_reject_later_source_change(tmp_path: Path, change: str):
    fixture = synthetic_source()
    snapshot = fixture.preflight(synthetic_config())
    original = tmp_path / "source.npy"
    fit = tmp_path / "fit.npy"
    values = np.vstack((fixture.fit.X.to_numpy(), fixture.test.X.to_numpy()))
    np.save(original, values)
    np.save(fit, fixture.fit.X.to_numpy())
    with original.open("rb") as handle:
        offset = inspect_npy_header_handle(handle, values.shape).data_offset
    digest = sha256_file(original)
    selected, stamp = _authenticated_source_rows(original, fit, snapshot, digest, offset)
    source = RealSource(tmp_path)
    source.config = {"inputs": {"all_gene_expression": "source.npy"}}
    source.full_provenance_verified = True
    source.authenticated_test_values = selected
    source.authenticated_source_stamp = stamp
    source._require_access = lambda _state: {"attempt_id": "synthetic"}
    snapshot = snapshot.__class__(
        snapshot.fit_pairs,
        snapshot.test_pairs,
        snapshot.gene_ids,
        snapshot.sample_ids,
        {**snapshot.input_hashes, "all_gene_expression": digest},
        snapshot.historical_config_compatibility,
    )
    if change == "in_place":
        changed = values.copy()
        changed[-1, -1] += 1
        np.save(original, changed)
    else:
        replacement = tmp_path / "replacement.npy"
        np.save(replacement, values)
        original.unlink()
        replacement.rename(original)
    with pytest.raises(ValueError, match="Authenticated source changed"):
        source._test_expression(snapshot)


@pytest.mark.parametrize("field", ["attempt_id", "authorization_sha256", "test_identity_sha256"])
def test_prediction_seal_binding_rejects_reuse_or_identity_change(tmp_path: Path, field: str):
    consume(tmp_path)
    seal_path = tmp_path / "result/prediction_seal.json"
    seal = json.loads(seal_path.read_text(encoding="utf-8"))
    seal[field] = "different-attempt-or-identity"
    seal_path.write_text(json.dumps(seal), encoding="utf-8")
    with pytest.raises(ValueError):
        verify_completed_test_only(tmp_path / "result")


def test_missing_terminal_evidence_and_failed_attempt_reuse_rejected(tmp_path: Path):
    source, config, path, fit_manifest = prepared(tmp_path)
    source.verify_full_provenance = lambda _snapshot: (_ for _ in ()).throw(
        ValueError("synthetic provenance failure")
    )
    with pytest.raises(ValueError, match="provenance failure"):
        execute_test_only(
            source,
            config,
            tmp_path / "failed",
            path,
            fit_manifest,
            authorization_sha256="synthetic",
            model_sha256=sha256_file(path),
            source_revision="synthetic",
        )
    with pytest.raises(FileExistsError):
        execute_test_only(
            source,
            config,
            tmp_path / "failed",
            path,
            fit_manifest,
            authorization_sha256="synthetic",
            model_sha256=sha256_file(path),
            source_revision="synthetic",
        )
    consume(tmp_path)
    (tmp_path / "result/success.json").unlink()
    with pytest.raises((FileNotFoundError, ValueError)):
        verify_completed_test_only(tmp_path / "result")


def test_production_rejects_caller_model_loader_before_access(tmp_path: Path):
    source, config, path, manifest = prepared(tmp_path)
    with pytest.raises(AuthorizationError, match="injected model loader"):
        execute_test_only(
            RealSource(tmp_path),
            config,
            tmp_path / "result",
            path,
            manifest,
            authorization_sha256="synthetic",
            model_sha256=sha256_file(path),
            source_revision="synthetic",
            attempt_id="synthetic",
            synthetic_model_loader=lambda _path: source,
        )
    assert not (tmp_path / "result").exists()


def test_direct_production_source_access_requires_authorization_and_claim(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    monkeypatch.setattr(consumer, "__file__", str(tmp_path / "src/oncorna/final_test_consumer.py"))
    source = RealSource(tmp_path)
    with pytest.raises(AuthorizationError, match="authorization, claim and output"):
        source.begin_test_access()
    with pytest.raises(AuthorizationError, match="authorization, claim and output"):
        source.load_test_labels(synthetic_source().preflight(synthetic_config()))


def test_completed_seal_detects_later_prediction_alteration(tmp_path: Path):
    consume(tmp_path)
    with (tmp_path / "result/test_predictions.tsv").open("a") as handle:
        handle.write("altered\n")
    with pytest.raises(ValueError):
        verify_completed_test_only(tmp_path / "result")


def test_label_access_before_prediction_seal_is_rejected(tmp_path: Path):
    fixture = synthetic_source()
    snapshot = fixture.preflight(synthetic_config())
    with pytest.raises(RuntimeError, match="before predictions"):
        fixture.load_test_labels(snapshot)
    with pytest.raises(RuntimeError, match="requires sealed predictions"):
        fixture.verify_label_provenance(snapshot)
    with pytest.raises(AuthorizationError, match="reviewed checkout"):
        RealSource(tmp_path).load_test_labels(snapshot)


@pytest.mark.parametrize("change", ["missing", "reordered", "duplicated", "gene_order"])
def test_identity_and_gene_order_fail_before_prediction(tmp_path: Path, change: str):
    source, config, path, fit_manifest = prepared(tmp_path)
    original = source.load_test_expression

    def altered(snapshot):
        test = original(snapshot)
        if change == "gene_order":
            return HeldOutData(
                test.X.loc[:, list(reversed(test.X.columns))], test.patient_ids, test.sample_ids
            )
        ids = list(test.sample_ids)
        patients = list(test.patient_ids)
        if change == "missing":
            ids = ids[:-1]
            patients = patients[:-1]
            frame = test.X.iloc[:-1]
        elif change == "reordered":
            ids.reverse()
            patients.reverse()
            frame = test.X.iloc[::-1]
        else:
            ids[1] = ids[0]
            patients[1] = patients[0]
            frame = test.X.copy()
            frame.index = ids
        return HeldOutData(frame, tuple(patients), tuple(ids))

    source.load_test_expression = altered
    with pytest.raises(ValueError):
        execute_test_only(
            source,
            config,
            tmp_path / "result",
            path,
            fit_manifest,
            authorization_sha256="synthetic",
            model_sha256=sha256_file(path),
            source_revision="synthetic",
        )
    assert "test_labels_loaded" not in source.events


def test_provenance_failure_stops_before_model_load_or_predictions(tmp_path: Path):
    source, config, path, fit_manifest = prepared(tmp_path)

    def failed(_snapshot):
        raise ValueError("Original-source checksum or row comparison failed")

    source.verify_full_provenance = failed

    def forbidden(_path):
        raise AssertionError("Model must not load after provenance failure")

    with pytest.raises(ValueError, match="Original-source"):
        execute_test_only(
            source,
            config,
            tmp_path / "result",
            path,
            fit_manifest,
            authorization_sha256="synthetic",
            model_sha256=sha256_file(path),
            source_revision="synthetic",
            synthetic_model_loader=forbidden,
        )
    assert "test_expression_loaded" not in source.events
    assert "test_labels_loaded" not in source.events


def test_changed_model_and_saved_feature_state_are_rejected(tmp_path: Path):
    source, config, path, fit_manifest = prepared(tmp_path)
    with pytest.raises(ValueError, match="Pipeline changed"):
        execute_test_only(
            source,
            config,
            tmp_path / "bad-hash",
            path,
            fit_manifest,
            authorization_sha256="synthetic",
            model_sha256="0" * 64,
            source_revision="synthetic",
        )
    assert "test_labels_loaded" not in source.events
    source2, config2, path2, manifest2 = prepared(tmp_path)
    manifest2["retained_gene_ids"] = list(reversed(manifest2["retained_gene_ids"]))
    with pytest.raises(ValueError, match="fitting state"):
        execute_test_only(
            source2,
            config2,
            tmp_path / "bad-features",
            path2,
            manifest2,
            authorization_sha256="synthetic",
            model_sha256=sha256_file(path2),
            source_revision="synthetic",
        )


def test_no_refitting_or_feature_selection_during_consumption(tmp_path: Path, monkeypatch):
    source, config, path, fit_manifest = prepared(tmp_path)
    pipeline = joblib.load(path)

    def forbidden(*_args, **_kwargs):
        raise AssertionError("A fitting method was called")

    for step in pipeline.named_steps.values():
        monkeypatch.setattr(step, "fit", forbidden)
        monkeypatch.setattr(step, "fit_transform", forbidden, raising=False)
    monkeypatch.setattr(pipeline, "fit", forbidden)
    execute_test_only(
        source,
        config,
        tmp_path / "result",
        path,
        fit_manifest,
        authorization_sha256="synthetic",
        model_sha256=sha256_file(path),
        source_revision="synthetic",
        synthetic_model_loader=lambda _path: pipeline,
    )


def test_seal_change_blocks_labels_and_postseal_change_blocks_completion(tmp_path: Path):
    source, config, path, fit_manifest = prepared(tmp_path)
    original = source.note_predictions_saved

    def tamper_before_labels(output, snapshot):
        with (tmp_path / "result/test_predictions.tsv").open("a") as handle:
            handle.write("altered\n")
        original(output, snapshot)

    source.note_predictions_saved = tamper_before_labels
    with pytest.raises(ValueError, match="Saved prediction count"):
        execute_test_only(
            source,
            config,
            tmp_path / "result",
            path,
            fit_manifest,
            authorization_sha256="synthetic",
            model_sha256=sha256_file(path),
            source_revision="synthetic",
        )
    assert "test_labels_loaded" not in source.events
    assert (tmp_path / "result/failure.json").is_file()


def test_authorization_contract_rejects_bad_attempt_or_bindings(tmp_path: Path):
    root = tmp_path.resolve()
    attempt = str(uuid.uuid4())
    record = expected_test_authorization(
        root, "a" * 40, "b" * 64, "c" * 64, attempt, "SYN-APPROVAL"
    )
    assert json.loads(_claim_payload(record, "d" * 64))["attempt_id"] == attempt
    with pytest.raises(ValueError):
        expected_test_authorization(
            root, "a" * 40, "b" * 64, "c" * 64, "wrong-attempt", "SYN-APPROVAL"
        )


def test_invalid_real_authorization_and_claim_rejected_without_access(
    tmp_path: Path,
    monkeypatch,
):
    root = tmp_path.resolve()
    (root / "configs").mkdir()
    (root / "data/processed/final_model_fit_v1").mkdir(parents=True)
    split = root / "data/processed/split_v1.json"
    split.write_text("{}\n")
    config = yaml.safe_load(
        (Path(__file__).resolve().parents[1] / "configs/final_test_v1.yaml").read_text()
    )
    config["inputs"]["expected_sha256"]["frozen_split"] = sha256_file(split)
    config_path = root / "configs/final_test_v1.yaml"
    config_path.write_text(yaml.safe_dump(config))
    fit_manifest = {"synthetic": True}
    (root / "data/processed/final_model_fit_v1/training_manifest.json").write_text(
        json.dumps(fit_manifest)
    )
    revision = "a" * 40

    def git_result(args, **_kwargs):
        return revision + "\n" if "rev-parse" in args else ""

    monkeypatch.setattr(consumer.subprocess, "check_output", git_result)
    monkeypatch.setattr(consumer, "verify_approved_fit", lambda _root: {})
    auth = root / consumer.AUTH_PATH
    record = expected_test_authorization(
        root,
        revision,
        sha256_file(config_path),
        sha256_file(split),
        str(uuid.uuid4()),
        "SYN-APPROVAL",
    )
    bad = copy.deepcopy(record)
    bad["fit_attempt_id"] = "incorrect-attempt"
    auth.write_text(json.dumps(bad))
    with pytest.raises(AuthorizationError, match="does not bind"):
        verify_test_authorization(root, auth)
    auth.write_text(json.dumps(record))
    accepted, auth_sha = verify_test_authorization(root, auth)
    assert accepted == record
    with pytest.raises(FileNotFoundError):
        consumer._verify_claimed_real_context(
            RealSource(root),
            root / consumer.OUTPUT_DIR,
            root / consumer.FIT_DIR / "pipeline.joblib",
            fit_manifest,
            config,
            auth_sha,
            consumer.MODEL_SHA256,
            revision,
        )
    claim = root / consumer.CLAIM_PATH
    claim.write_text("wrong claim")

    def reject_claim(_path, _payload):
        raise AuthorizationError("Invalid test attempt claim")

    monkeypatch.setattr(consumer, "_verify_existing_claim", reject_claim)
    with pytest.raises(AuthorizationError, match="Invalid test attempt claim"):
        consumer._verify_claimed_real_context(
            RealSource(root),
            root / consumer.OUTPUT_DIR,
            root / consumer.FIT_DIR / "pipeline.joblib",
            fit_manifest,
            config,
            auth_sha,
            consumer.MODEL_SHA256,
            revision,
        )
    assert not (root / consumer.OUTPUT_DIR).exists()


@pytest.mark.parametrize("changed", ["model", "manifest"])
def test_offline_fit_identities_cannot_be_replaced_by_manifest(
    tmp_path: Path,
    monkeypatch,
    changed: str,
):
    output = tmp_path / "data/processed/final_model_fit_v1"
    output.mkdir(parents=True)
    (output / "pipeline.joblib").write_bytes(b"other model")
    (output / "manifest.json").write_text("{}")
    if changed == "manifest":
        monkeypatch.setattr(consumer, "MODEL_SHA256", sha256_file(output / "pipeline.joblib"))
    with pytest.raises(AuthorizationError, match="offline identities"):
        verify_approved_fit(tmp_path)


@pytest.mark.parametrize("changed", ["attempt_id", "authorization_sha256"])
def test_incorrect_fitting_attempt_or_authorization_rejected(
    tmp_path: Path,
    monkeypatch,
    changed: str,
):
    output = tmp_path / consumer.FIT_DIR
    output.mkdir(parents=True)
    (output / "pipeline.joblib").write_bytes(b"synthetic pipeline")
    (output / "manifest.json").write_text("{}")
    monkeypatch.setattr(consumer, "MODEL_SHA256", sha256_file(output / "pipeline.joblib"))
    monkeypatch.setattr(consumer, "MANIFEST_SHA256", sha256_file(output / "manifest.json"))
    manifest = {
        "scope": "production",
        "source_revision": consumer.FIT_REVISION,
        "attempt_id": consumer.FIT_ATTEMPT,
        "authorization_sha256": consumer.FIT_AUTH_SHA256,
        "output_hashes": {"pipeline.joblib": consumer.MODEL_SHA256},
        "fit_count": 675,
        "original_gene_count": 20530,
        "retained_gene_count": 17623,
    }
    manifest[changed] = "incorrect"
    monkeypatch.setattr(consumer, "verify_sealed_fit", lambda *_args, **_kwargs: manifest)
    with pytest.raises(AuthorizationError, match="approved production attempt"):
        verify_approved_fit(tmp_path)
