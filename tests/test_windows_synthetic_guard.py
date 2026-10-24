"""Windows must reject synthetic filesystem execution before touching candidate paths."""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import oncorna.final_test as legacy
import oncorna.final_test_consumer as consumer
import oncorna.synthetic_paths as paths
from oncorna.final_test import FitData, SyntheticFileSource, SyntheticFixtureSource
from oncorna.final_test import TestData as HeldOutData
from oncorna.synthetic_paths import SyntheticIsolationError

pytestmark = pytest.mark.skipif(os.name != "nt", reason="Windows production checkout boundary")


def fixture_source() -> SyntheticFixtureSource:
    fit = FitData(pd.DataFrame(), np.asarray([], dtype=object), (), ())
    test = HeldOutData(pd.DataFrame(), (), ())
    return SyntheticFixtureSource(fit, test, np.asarray([], dtype=object))


def forbidden(*_args: object, **_kwargs: object) -> None:
    pytest.fail("Synthetic filesystem side effect reached")


def test_consumer_rejects_before_output_input_or_model_load(tmp_path: Path, monkeypatch):
    source = fixture_source()
    output = tmp_path / "result"
    model = tmp_path / "missing-model.joblib"
    monkeypatch.setattr(consumer, "RunState", forbidden)
    monkeypatch.setattr(consumer, "synthetic_output_path", forbidden)
    monkeypatch.setattr(consumer, "synthetic_input_path", forbidden)
    monkeypatch.setattr(consumer, "sha256_file", forbidden)
    monkeypatch.setattr(consumer.joblib, "load", forbidden)
    monkeypatch.setattr(source, "preflight", forbidden)
    with pytest.raises(SyntheticIsolationError):
        consumer.execute_test_only(
            source,
            {},
            output,
            model,
            {},
            authorization_sha256="synthetic",
            model_sha256="synthetic",
            source_revision="synthetic",
            synthetic_model_loader=forbidden,
        )
    assert not output.exists()
    assert not model.exists()


def test_legacy_executor_rejects_before_preflight_or_output(tmp_path: Path, monkeypatch):
    source = fixture_source()
    output = tmp_path / "legacy-result"
    monkeypatch.setattr(legacy, "RunState", forbidden)
    monkeypatch.setattr(legacy, "validate_protocol", forbidden)
    monkeypatch.setattr(source, "preflight", forbidden)
    with pytest.raises(SyntheticIsolationError):
        legacy.execute_once(
            source,
            {},
            output,
            config_sha256="synthetic",
            source_revision="synthetic",
            authorization_sha256=None,
        )
    assert not output.exists()


def test_file_source_constructor_rejects_before_path_inspection(tmp_path: Path, monkeypatch):
    source = fixture_source()
    labels = tmp_path / "missing-labels.tsv"
    monkeypatch.setattr(legacy, "synthetic_input_path", forbidden)
    with pytest.raises(SyntheticIsolationError):
        SyntheticFileSource(source.fit, source.test, labels)
    assert not labels.exists()


def test_direct_path_helpers_reject_before_file_operations(tmp_path: Path, monkeypatch):
    candidate = tmp_path / "missing-input"
    with monkeypatch.context() as context:
        context.setattr(Path, "resolve", forbidden)
        context.setattr(Path, "stat", forbidden)
        context.setattr(Path, "open", forbidden)
        for operation in (
            lambda: paths.synthetic_output_path(candidate),
            lambda: paths.synthetic_input_path(candidate),
            lambda: paths.verify_synthetic_input(candidate, (0, 0, 0, 0, 0)),
            lambda: paths.sha256_synthetic_input(candidate, (0, 0, 0, 0, 0)),
        ):
            with pytest.raises(SyntheticIsolationError):
                operation()
    assert not candidate.exists()


def test_real_entrypoint_still_requires_separate_approval(tmp_path: Path, monkeypatch):
    monkeypatch.setattr(consumer, "verify_approved_fit", lambda _root: {})

    def reject_authorization(*_args: object) -> None:
        raise legacy.AuthorizationError("Separate test authorization required")

    monkeypatch.setattr(consumer, "verify_test_authorization", reject_authorization)
    monkeypatch.setattr(consumer, "RealSource", forbidden)
    with pytest.raises(legacy.AuthorizationError, match="Separate test authorization"):
        consumer.run_real(tmp_path, tmp_path / consumer.AUTH_PATH)
    assert not (tmp_path / consumer.OUTPUT_DIR).exists()
    assert not (tmp_path / consumer.CLAIM_PATH).exists()
