"""Synthetic checks for test-only model consumption and irreversible ordering."""

from __future__ import annotations

import copy
import json
import os
import uuid
from pathlib import Path

import joblib
import numpy as np
import pytest
import yaml
from test_final_test import synthetic_config, synthetic_source

import oncorna.final_test as final_test
import oncorna.final_test_consumer as consumer
from oncorna.final_fit_export import inspect_npy_header_handle
from oncorna.final_test import (
    LABELS,
    AuthorizationError,
    RealSource,
    Snapshot,
    SyntheticFileSource,
    _authenticated_source_rows,
    _write_predictions,
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
from oncorna.preprocessing import identity_digest, sha256_file
from oncorna.synthetic_paths import (
    PROJECT_DATA,
    synthetic_input_path,
    synthetic_output_path,
    verify_synthetic_input,
)

pytestmark = pytest.mark.skipif(
    os.name == "nt",
    reason="Filesystem-backed synthetic integration runs in the code-only container",
)


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


@pytest.mark.parametrize(
    "reserved_name",
    [
        "final_test_v1",
        "final_test_v1_authorization.json",
        "final_test_v1_attempt.json",
    ],
)
def test_synthetic_cannot_reserve_real_test_output(tmp_path: Path, monkeypatch, reserved_name: str):
    source, config, model, fit_manifest = prepared(tmp_path)
    reserved = PROJECT_DATA / "processed" / reserved_name
    assert not reserved.exists()
    monkeypatch.setattr(
        consumer, "RunState", lambda *_args, **_kwargs: pytest.fail("mkdir reached")
    )
    with pytest.raises(ValueError, match="protected production data"):
        execute_test_only(
            source,
            config,
            reserved,
            model,
            fit_manifest,
            authorization_sha256="synthetic",
            model_sha256=sha256_file(model),
            source_revision="synthetic",
        )
    assert not reserved.exists()


def test_synthetic_manifest_and_seal_have_distinct_scope(tmp_path: Path):
    manifest, _ = consume(tmp_path)
    seal = json.loads((tmp_path / "result/prediction_seal.json").read_text())
    assert manifest["scope"] == seal["scope"] == "synthetic_test_only"
    assert verify_completed_test_only(tmp_path / "result")["scope"] == "synthetic_test_only"


def test_synthetic_paths_reject_case_variants_before_open(tmp_path: Path, monkeypatch):
    protected = tmp_path / "data"
    location = protected / "processed"
    location.mkdir(parents=True)
    candidate = location / "synthetic_labels.tsv"
    candidate.write_text("synthetic", encoding="utf-8")
    upper = tmp_path / "DATA/PROCESSED/synthetic_labels.tsv"
    if os.name == "nt":
        assert upper.samefile(candidate)
    opened = []
    original_open = Path.open

    def watched_open(path, *args, **kwargs):
        if path.name == candidate.name:
            opened.append(path)
            pytest.fail("Protected candidate opened")
        return original_open(path, *args, **kwargs)

    monkeypatch.setattr(Path, "open", watched_open)
    for path in (candidate, upper):
        with pytest.raises(ValueError):
            synthetic_input_path(path, protected_data=protected)
    assert opened == []


def test_file_source_rejects_protected_label_before_open(tmp_path: Path, monkeypatch):
    fixture = synthetic_source()
    protected = tmp_path / "data"
    location = protected / "processed"
    location.mkdir(parents=True)
    candidate = location / "labels.tsv"
    candidate.write_text("synthetic", encoding="utf-8")
    original_check = final_test.synthetic_input_path
    monkeypatch.setattr(
        final_test,
        "synthetic_input_path",
        lambda path: original_check(path, protected_data=protected),
    )
    monkeypatch.setattr(Path, "open", lambda *_args, **_kwargs: pytest.fail("opened"))
    with pytest.raises(ValueError, match="protected production data"):
        SyntheticFileSource(fixture.fit, fixture.test, candidate)


def test_synthetic_alias_and_changed_path_are_rejected(tmp_path: Path):
    protected = tmp_path / "data"
    protected.mkdir()
    target = protected / "protected.tsv"
    target.write_text("synthetic", encoding="utf-8")
    alias = tmp_path / "alias.tsv"
    directory_alias = tmp_path / "data_alias"
    try:
        alias.symlink_to(target)
        directory_alias.symlink_to(protected, target_is_directory=True)
    except (OSError, NotImplementedError):
        pytest.skip("File symlinks require privileges on this Windows host")
    with pytest.raises(ValueError, match="symlink|protected"):
        synthetic_input_path(alias, protected_data=protected)
    with pytest.raises(ValueError, match="symlink|protected"):
        synthetic_output_path(directory_alias / "child", protected_data=protected)


def test_synthetic_path_substitution_rejected(tmp_path: Path):
    protected = tmp_path / "data"
    protected.mkdir()
    allowed = tmp_path / "allowed.tsv"
    allowed.write_text("first", encoding="utf-8")
    stamp = synthetic_input_path(allowed, protected_data=protected)
    allowed.unlink()
    allowed.write_text("second", encoding="utf-8")
    with pytest.raises(ValueError, match="changed"):
        verify_synthetic_input(allowed, stamp, protected_data=protected)


def test_label_substitution_is_rejected_before_open(tmp_path: Path, monkeypatch):
    fixture, config, model, fit_manifest = prepared(tmp_path)
    labels = tmp_path / "labels.tsv"
    labels.write_text("patient_id\tsample_id\ttrue_label\n", encoding="utf-8")
    source = SyntheticFileSource(fixture.fit, fixture.test, labels)
    labels.unlink()
    labels.write_text("replacement", encoding="utf-8")
    original_open = Path.open

    def watched_open(path, *args, **kwargs):
        if path == labels:
            pytest.fail("Changed label path opened")
        return original_open(path, *args, **kwargs)

    monkeypatch.setattr(Path, "open", watched_open)
    with pytest.raises(ValueError, match="changed"):
        execute_test_only(
            source,
            config,
            tmp_path / "result",
            model,
            fit_manifest,
            authorization_sha256="synthetic",
            model_sha256=sha256_file(model),
            source_revision="synthetic",
        )
    assert "label_file_opened" not in source.events


def test_synthetic_model_path_checked_before_hash_or_load(tmp_path: Path, monkeypatch):
    source, config, model, fit_manifest = prepared(tmp_path)
    protected = tmp_path / "data"
    protected.mkdir()
    candidate = protected / "synthetic.joblib"
    model.replace(candidate)
    import oncorna.synthetic_paths as paths

    original_check = paths.synthetic_input_path
    monkeypatch.setattr(
        consumer,
        "synthetic_input_path",
        lambda path: original_check(path, protected_data=protected),
    )
    monkeypatch.setattr(
        consumer, "RunState", lambda *_args, **_kwargs: pytest.fail("mkdir reached")
    )
    with pytest.raises(ValueError, match="protected production data"):
        execute_test_only(
            source,
            config,
            tmp_path / "result",
            candidate,
            fit_manifest,
            authorization_sha256="synthetic",
            model_sha256="0" * 64,
            source_revision="synthetic",
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


def production_seal_fixture(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    """Use the real source seal verifier with only disposable synthetic files."""
    root = tmp_path.resolve()
    output = root / consumer.OUTPUT_DIR
    output.mkdir(parents=True)
    cohort = root / "synthetic_cohort.tsv"
    cohort.write_text("patient_id\tsample_id\tpam50_normalized_label\nSYN-P1\tSYN-S1\tLuminal A\n")
    source = RealSource(root)
    source.config = {"inputs": {"cohort": cohort.name}}
    source.full_provenance_verified = True
    snapshot = Snapshot((), (("SYN-P1", "SYN-S1"),), (), (), {})
    approved = expected_test_authorization(
        root, "a" * 40, "b" * 64, "c" * 64, str(uuid.uuid4()), "SYN-APPROVAL"
    )
    approved["authorization_sha256"] = "d" * 64

    def synthetic_approval(checked_root: Path, required_state: str):
        assert checked_root == root
        assert required_state == "predictions_sealed"
        return approved

    monkeypatch.setattr(consumer, "verify_active_test_attempt", synthetic_approval)
    rows = [
        {
            "patient_id": "SYN-P1",
            "sample_id": "SYN-S1",
            "predicted_label": "Luminal A",
            **{
                f"probability_{label}": value
                for label, value in zip(LABELS, (0.7, 0.3, 0.0, 0.0, 0.0), strict=True)
            },
        }
    ]
    prediction_path = output / "test_predictions.tsv"
    _write_predictions(prediction_path, rows, LABELS)
    identity_sha = identity_digest(snapshot.test_pairs)
    seal = {
        "status": "sealed_before_labels",
        "scope": "one_shot_final_test",
        "attempt_id": approved["attempt_id"],
        "authorization_sha256": approved["authorization_sha256"],
        "source_revision": approved["approved_commit"],
        "model_sha256": approved["pipeline_sha256"],
        "test_identity_sha256": identity_sha,
        "prediction_sha256": sha256_file(prediction_path),
        "row_count": 1,
        "class_order": list(LABELS),
    }
    seal_path = output / "prediction_seal.json"
    seal_path.write_text(json.dumps(seal))
    entries = [
        {"sequence": 1, "state": "reserved"},
        {"sequence": 2, "state": "test_access_started"},
        {
            "sequence": 3,
            "state": "predictions_sealed",
            "prediction_sha256": seal["prediction_sha256"],
            "seal_sha256": sha256_file(seal_path),
            "test_identity_sha256": identity_sha,
        },
    ]
    state_path = output / "run_state.jsonl"
    state_path.write_text("".join(json.dumps(entry) + "\n" for entry in entries))
    return source, snapshot, approved, output, cohort, seal, entries


def test_production_seal_fixture_passes_before_synthetic_label_open(tmp_path: Path, monkeypatch):
    source, snapshot, _approved, output, cohort, _seal, _entries = production_seal_fixture(
        tmp_path, monkeypatch
    )
    opened = []
    original_open = Path.open

    def watch_label_open(path, *args, **kwargs):
        if path == cohort:
            opened.append(path)
        return original_open(path, *args, **kwargs)

    monkeypatch.setattr(Path, "open", watch_label_open)
    with pytest.raises(AuthorizationError, match="sealed predictions"):
        source.load_test_labels(snapshot)
    assert opened == []
    source.note_predictions_saved(output, snapshot)
    assert source.predictions_saved
    assert source.load_test_labels(snapshot).tolist() == ["Luminal A"]
    assert opened == [cohort]
    assert not (tmp_path / consumer.AUTH_PATH).exists()
    assert not (tmp_path / consumer.CLAIM_PATH).exists()
    with pytest.raises(ValueError, match="terminal completed state"):
        verify_completed_test_only(output)


@pytest.mark.parametrize(
    "change",
    [
        "approved_revision",
        "pipeline_digest",
        "attempt_id",
        "authorization_sha256",
        "missing_scope",
        "wrong_scope",
        "seal_attempt_id",
        "seal_authorization_sha256",
        "test_identity_digest",
        "class_order",
        "unexpected_seal_field",
        "missing_prediction_sha256",
        "predicted_label",
        "probability",
        "audit_state",
    ],
)
def test_production_seal_fixture_rejects_before_label_open(tmp_path: Path, monkeypatch, change):
    source, snapshot, approved, output, cohort, seal, entries = production_seal_fixture(
        tmp_path, monkeypatch
    )
    seal_path = output / "prediction_seal.json"
    prediction_path = output / "test_predictions.tsv"
    if change == "approved_revision":
        approved["approved_commit"] = "e" * 40
    elif change == "pipeline_digest":
        approved["pipeline_sha256"] = "f" * 64
    elif change == "attempt_id":
        approved["attempt_id"] = str(uuid.uuid4())
    elif change == "authorization_sha256":
        approved["authorization_sha256"] = "0" * 64
    elif change in (
        "missing_scope",
        "wrong_scope",
        "seal_attempt_id",
        "seal_authorization_sha256",
        "test_identity_digest",
        "class_order",
        "unexpected_seal_field",
        "missing_prediction_sha256",
    ):
        if change == "missing_scope":
            del seal["scope"]
        elif change == "wrong_scope":
            seal["scope"] = "synthetic_test_only"
        elif change == "seal_attempt_id":
            seal["attempt_id"] = str(uuid.uuid4())
        elif change == "seal_authorization_sha256":
            seal["authorization_sha256"] = "0" * 64
        elif change == "test_identity_digest":
            seal["test_identity_sha256"] = "0" * 64
        elif change == "class_order":
            seal["class_order"] = list(reversed(LABELS))
        elif change == "unexpected_seal_field":
            seal["unexpected"] = "SYN-EXTRA"
        else:
            del seal["prediction_sha256"]
        seal_path.write_text(json.dumps(seal))
        entries[-1]["seal_sha256"] = sha256_file(seal_path)
    elif change in ("predicted_label", "probability"):
        original = prediction_path.read_text()
        altered = (
            original.replace("SYN-S1\tLuminal A\t", "SYN-S1\tLuminal B\t")
            if change == "predicted_label"
            else original.replace("\t0.7\t0.3\t", "\t0.6\t0.4\t")
        )
        assert altered != original
        prediction_path.write_text(altered)
    else:
        entries[-1]["state"] = "completed"
    (output / "run_state.jsonl").write_text("".join(json.dumps(entry) + "\n" for entry in entries))
    opened = []
    original_open = Path.open

    def watch_label_open(path, *args, **kwargs):
        if path == cohort:
            opened.append(path)
        return original_open(path, *args, **kwargs)

    monkeypatch.setattr(Path, "open", watch_label_open)
    with pytest.raises(AuthorizationError, match="Prediction seal differs"):
        source.note_predictions_saved(output, snapshot)
    assert not source.predictions_saved
    with pytest.raises(AuthorizationError):
        source.load_test_labels(snapshot)
    assert opened == []


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
