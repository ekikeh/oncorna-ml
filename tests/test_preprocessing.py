from __future__ import annotations

import gzip
import hashlib
import json

import numpy as np
import pandas as pd
import pytest
from sklearn.base import clone
from sklearn.model_selection import KFold
from sklearn.pipeline import Pipeline

from oncorna.preprocessing import (
    AllGeneView,
    GenePrevalenceFilter,
    partition_frames,
    read_all_gene_view,
    sha256_file,
    validate_frozen_partitions,
    verify_historical_config,
)


def _frame(values: list[list[float]], genes: list[str]) -> pd.DataFrame:
    return pd.DataFrame(
        values,
        index=[f"S{index}" for index in range(len(values))],
        columns=genes,
    )


def test_strict_threshold_ceiling_and_gene_order() -> None:
    X = _frame(
        [
            [1.0, 1.1, 2.0],
            [1.1, np.nan, 1.1],
            [2.0, 1.1, 1.0],
            [np.nan, 0.0, 1.1],
            [0.0, 0.0, 0.0],
        ],
        ["boundary", "too_rare", "retained"],
    )
    fitted = GenePrevalenceFilter(expression_gt=1.0, minimum_sample_fraction=0.41).fit(X)
    assert fitted.required_count_ == 3
    assert fitted.prevalence_counts_.tolist() == [2, 2, 3]
    assert fitted.selected_gene_ids_ == ("retained",)
    assert fitted.get_support().tolist() == [False, False, True]
    assert fitted.get_feature_names_out().tolist() == ["retained"]
    assert fitted.missing_value_count_ == 2
    assert fitted.transform(X).index.tolist() == X.index.tolist()


def test_full_training_count_is_102_and_labels_do_not_affect_mask() -> None:
    values = np.zeros((506, 2), dtype=np.float64)
    values[:101, 0] = 2
    values[:102, 1] = 2
    X = pd.DataFrame(values, index=[f"S{n}" for n in range(506)], columns=["low", "pass"])
    first = GenePrevalenceFilter().fit(X, y=np.zeros(506))
    second = GenePrevalenceFilter().fit(X, y=np.arange(506) % 5)
    assert first.required_count_ == 102
    assert first.selected_gene_ids_ == second.selected_gene_ids_ == ("pass",)
    assert first.fitting_sample_ids_ == tuple(X.index)


def test_infinite_values_are_rejected_on_fit_and_transform() -> None:
    X = _frame([[2.0, np.nan], [2.0, 3.0]], ["A", "B"])
    fitted = GenePrevalenceFilter(minimum_sample_fraction=0.5).fit(X)
    assert fitted.missing_value_count_ == 1
    changed = X.copy()
    changed.loc["S0", "A"] = np.inf
    with pytest.raises(ValueError, match="Infinite expression"):
        GenePrevalenceFilter().fit(changed)
    with pytest.raises(ValueError, match="Infinite expression"):
        fitted.transform(changed)


def test_transform_reuses_mask_and_rejects_schema_changes() -> None:
    train = _frame([[2.0, 0.0], [2.0, 0.0], [0.0, 0.0]], ["A", "B"])
    held_out = _frame([[0.0, 500.0], [0.0, 500.0]], ["A", "B"])
    held_out.index = ["V0", "T0"]
    fitted = GenePrevalenceFilter(minimum_sample_fraction=0.5).fit(train)
    original_mask = fitted.get_support()
    transformed = fitted.transform(held_out)
    assert transformed.columns.tolist() == ["A"]
    assert transformed.index.tolist() == held_out.index.tolist()
    assert np.array_equal(fitted.get_support(), original_mask)
    with pytest.raises(ValueError, match="schema"):
        fitted.transform(held_out[["B", "A"]])
    with pytest.raises(ValueError, match="schema"):
        fitted.transform(held_out.rename(columns={"A": "different"}))
    with pytest.raises(ValueError, match="feature names"):
        fitted.get_feature_names_out(["B", "A"])


def test_repeated_fits_are_deterministic_and_input_must_have_gene_ids() -> None:
    X = _frame([[2.0, 0.0], [3.0, 2.0]], ["G2", "G1"])
    first = GenePrevalenceFilter().fit(X)
    second = GenePrevalenceFilter().fit(X)
    assert first.selected_gene_ids_ == second.selected_gene_ids_ == ("G2", "G1")
    assert np.array_equal(first.prevalence_counts_, second.prevalence_counts_)
    with pytest.raises(TypeError, match="DataFrame"):
        GenePrevalenceFilter().fit(X.to_numpy())  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="unique"):
        GenePrevalenceFilter().fit(X.set_axis(["G2", "G2"], axis=1))


def test_pipeline_refits_filter_in_each_synthetic_cv_fold() -> None:
    X = _frame(
        [
            [2.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [0.0, 2.0, 0.0],
            [0.0, 2.0, 0.0],
            [0.0, 0.0, 2.0],
            [0.0, 0.0, 2.0],
        ],
        ["A", "B", "C"],
    )
    pipeline = Pipeline([("gene_filter", GenePrevalenceFilter(minimum_sample_fraction=0.25))])
    masks: list[tuple[str, ...]] = []
    for fit_indices, held_out_indices in KFold(n_splits=3).split(X):
        fitted = clone(pipeline).fit(X.iloc[fit_indices])
        gene_filter = fitted.named_steps["gene_filter"]
        assert set(gene_filter.fitting_sample_ids_).isdisjoint(X.iloc[held_out_indices].index)
        assert gene_filter.fitting_sample_ids_ == tuple(X.iloc[fit_indices].index)
        assert tuple(fitted.transform(X.iloc[held_out_indices]).columns) == (
            gene_filter.selected_gene_ids_
        )
        masks.append(gene_filter.selected_gene_ids_)
    assert masks == [("B", "C"), ("A", "C"), ("A", "B")]


def test_all_gene_reader_and_partition_frames_preserve_order(tmp_path) -> None:
    source = tmp_path / "source.tsv.gz"
    with gzip.open(source, "wt", encoding="utf-8", newline="") as handle:
        handle.write("sample\tS3\tS1\tS2\n")
        handle.write("G2\t9\t1\t2\n")
        handle.write("G1\t8\tNA\t3\n")
    view = read_all_gene_view(source, ["S1", "S2", "S3"], expected_gene_count=2)
    assert view.gene_ids == ("G2", "G1")
    assert view.sample_ids == ("S1", "S2", "S3")
    assert view.source_sample_count == 3
    assert view.values.dtype == np.float64
    assert view.values[0, 0] == 1.0
    assert np.isnan(view.values[0, 1])
    parts = {
        "train": (("P2", "S2"),),
        "validation": (("P1", "S1"),),
        "test": (("P3", "S3"),),
    }
    frames = partition_frames(view, parts)
    assert frames["train"].columns.tolist() == ["G2", "G1"]
    assert frames["train"].index.tolist() == ["S2"]
    assert frames["train"].iloc[0].tolist() == [2.0, 3.0]


def test_only_frozen_training_patient_sample_pairs_enter_fit() -> None:
    records, split = _synthetic_split()
    partitions = validate_frozen_partitions(records, split)
    view = AllGeneView(
        gene_ids=("A", "B"),
        sample_ids=("S0", "S1", "S2"),
        values=np.asarray([[2.0, 0.0], [0.0, 9.0], [0.0, 9.0]]),
        source_sample_count=3,
    )
    frames = partition_frames(view, partitions)
    fitted = GenePrevalenceFilter().fit(frames["train"])
    held_out_patients = {
        patient for name in ("validation", "test") for patient, _ in partitions[name]
    }
    training_patients = {patient for patient, _ in partitions["train"]}
    assert training_patients.isdisjoint(held_out_patients)
    assert fitted.fitting_sample_ids_ == ("S0",)
    assert fitted.selected_gene_ids_ == ("A",)
    assert fitted.transform(frames["validation"]).columns.tolist() == ["A"]
    assert fitted.transform(frames["test"]).columns.tolist() == ["A"]


def test_all_gene_reader_rejects_infinity(tmp_path) -> None:
    source = tmp_path / "source.tsv.gz"
    with gzip.open(source, "wt", encoding="utf-8") as handle:
        handle.write("sample\tS1\nG1\tinf\n")
    with pytest.raises(ValueError, match="Infinite expression"):
        read_all_gene_view(source, ["S1"], expected_gene_count=1)


def _synthetic_split() -> tuple[list[dict[str, str]], dict]:
    records = [
        {"patient_id": f"P{i}", "sample_id": f"S{i}", "pam50_normalized_label": "A"}
        for i in range(3)
    ]
    parts = {}
    for name, i in (("train", 0), ("validation", 1), ("test", 2)):
        parts[name] = {
            "patient_ids": [f"P{i}"],
            "sample_ids": [f"S{i}"],
            "patient_sample_map": [{"patient_id": f"P{i}", "sample_id": f"S{i}"}],
            "n_patients": 1,
            "n_samples": 1,
            "pam50_counts": {"A": 1},
        }
    split = {
        "schema_version": 1,
        "split_id": "split_v1",
        "split_unit": "patient_id",
        "validation_results": {"passed": True},
        "source_cohort": {"patient_count": 3, "sample_count": 3},
        "cohort_pam50_counts": {"A": 3},
        "partitions": parts,
    }
    return records, split


def test_frozen_split_validation_is_read_only_and_rejects_overlap(tmp_path) -> None:
    records, split = _synthetic_split()
    path = tmp_path / "split_v1.json"
    path.write_text(json.dumps(split), encoding="utf-8")
    before = sha256_file(path)
    result = validate_frozen_partitions(records, json.loads(path.read_text(encoding="utf-8")))
    assert result["train"] == (("P0", "S0"),)
    assert sha256_file(path) == before
    assert path.read_bytes() == json.dumps(split).encode("utf-8")
    split["partitions"]["test"] = split["partitions"]["train"]
    with pytest.raises(ValueError, match="overlap"):
        validate_frozen_partitions(records, split)


def test_config_line_endings_are_checked_without_changing_recorded_hash(tmp_path) -> None:
    config = tmp_path / "default.yaml"
    lf = b"random_seed: 20261006\nsplitting:\n  split_unit: patient_id\n"
    config.write_bytes(lf)
    recorded = hashlib.sha256(lf).hexdigest()
    assert verify_historical_config(config, recorded)["compatibility_mode"] == "exact_raw_bytes"
    config.write_bytes(lf.replace(b"\n", b"\r\n"))
    result = verify_historical_config(config, recorded)
    assert result["compatibility_mode"] == "recorded_lf_current_crlf"
    assert result["recorded_raw_sha256"] == recorded
    config.write_bytes(b"random_seed: 1\r\n")
    with pytest.raises(ValueError, match="beyond LF/CRLF"):
        verify_historical_config(config, recorded)
