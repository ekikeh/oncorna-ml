"""Synthetic-only coverage of real fitting metadata and authorization boundaries."""

from __future__ import annotations

import json
import uuid
from pathlib import Path

import pytest
import yaml
from test_final_test import synthetic_real_adapter_metadata

import oncorna.final_model_fit_real as real
from oncorna.final_fit_export import sequence_digest
from oncorna.final_model_fit import (
    ATTEMPT_PATH,
    AUTHORIZATION_PATH,
    OUTPUT_DIRECTORY,
    _digest,
    expected_authorization,
    run_real_fit,
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
