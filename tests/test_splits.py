from __future__ import annotations

import csv
import gzip
import hashlib
import json
from collections import Counter
from datetime import datetime
from pathlib import Path

import pytest

from oncorna.cohort import PAM50_CLASS_ORDER
from oncorna.splits import (
    generate_split_manifest,
    largest_remainder_sizes,
    write_frozen_split_manifest,
)

COHORT_FIELDS = (
    "sample_id",
    "patient_id",
    "pam50_original_label",
    "pam50_normalized_label",
    "sample_type",
)
PROPORTIONS = {"train": 0.60, "validation": 0.20, "test": 0.20}
TIE_ORDER = ("validation", "test", "train")


def _write_cohort(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=COHORT_FIELDS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def _write_matrix(path: Path, sample_ids: list[str], values: str = "not-a-number") -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["gene_id", *sample_ids])
        # Deliberately nonnumeric: split generation must never parse expression rows.
        writer.writerow(["g1", *([values] * len(sample_ids))])


@pytest.fixture
def split_project(tmp_path: Path) -> dict[str, Path | list[dict[str, str]]]:
    records: list[dict[str, str]] = []
    for class_index, label in enumerate(PAM50_CLASS_ORDER):
        for sample_index in range(10):
            index = class_index * 10 + sample_index
            records.append(
                {
                    "sample_id": f"S{index:03d}",
                    "patient_id": f"P{index:03d}",
                    "pam50_original_label": label,
                    "pam50_normalized_label": label,
                    "sample_type": "Primary Tumor",
                }
            )
    cohort_path = tmp_path / "data/processed/classification_cohort.tsv"
    _write_cohort(cohort_path, records)
    sample_ids = [row["sample_id"] for row in records]

    matrix_path = tmp_path / "data/processed/ml_expression_matrix.tsv.gz"
    _write_matrix(matrix_path, sample_ids)
    matrix_manifest_path = tmp_path / "data/processed/ml_matrix_manifest.json"
    matrix_manifest_path.write_text(
        json.dumps(
            {
                "expression_units": "log2(normalized_count + 1)",
                "sample_order": sample_ids,
                "matrix": {
                    "path": "data/processed/ml_expression_matrix.tsv.gz",
                    "sha256": "0" * 64,
                    "shape": {
                        "gene_rows": 1,
                        "sample_columns": len(sample_ids),
                        "total_columns_including_gene_id": len(sample_ids) + 1,
                    },
                },
            }
        ),
        encoding="utf-8",
    )
    config_path = tmp_path / "configs/default.yaml"
    config_path.parent.mkdir(parents=True, exist_ok=True)
    config_path.write_text("random_seed: 20261006\n", encoding="utf-8")
    return {
        "root": tmp_path,
        "cohort": cohort_path,
        "matrix": matrix_path,
        "matrix_manifest": matrix_manifest_path,
        "config": config_path,
        "output": tmp_path / "data/processed/split_v1.json",
        "records": records,
    }


def _generate(
    paths: dict[str, Path | list[dict[str, str]]],
    *,
    seed: int = 20261006,
    proportions: dict[str, float] | None = None,
) -> dict:
    return generate_split_manifest(
        cohort_path=paths["cohort"],  # type: ignore[arg-type]
        expression_matrix_path=paths["matrix"],  # type: ignore[arg-type]
        matrix_manifest_path=paths["matrix_manifest"],  # type: ignore[arg-type]
        config_path=paths["config"],  # type: ignore[arg-type]
        project_root=paths["root"],  # type: ignore[arg-type]
        expected_patient_count=50,
        expected_gene_count=1,
        expression_scale="log2(normalized_count + 1)",
        random_seed=seed,
        proportions=proportions or PROPORTIONS,
        rounding_tie_order=TIE_ORDER,
        test_random_state_offset=0,
        validation_random_state_offset=1,
    )


def test_patient_and_sample_partitions_are_disjoint_and_complete(
    split_project: dict[str, Path | list[dict[str, str]]],
) -> None:
    manifest = _generate(split_project)
    records = split_project["records"]
    assert isinstance(records, list)
    expected_patients = {row["patient_id"] for row in records}
    expected_samples = {row["sample_id"] for row in records}

    patient_sets = {
        split: set(manifest["partitions"][split]["patient_ids"])
        for split in ("train", "validation", "test")
    }
    sample_sets = {
        split: set(manifest["partitions"][split]["sample_ids"])
        for split in ("train", "validation", "test")
    }
    for left, right in (("train", "validation"), ("train", "test"), ("validation", "test")):
        assert patient_sets[left].isdisjoint(patient_sets[right])
        assert sample_sets[left].isdisjoint(sample_sets[right])
    assert set.union(*patient_sets.values()) == expected_patients
    assert set.union(*sample_sets.values()) == expected_samples

    patient_to_sample = {row["patient_id"]: row["sample_id"] for row in records}
    observed_map: dict[str, str] = {}
    for split in ("train", "validation", "test"):
        partition = manifest["partitions"][split]
        assert partition["n_patients"] == partition["n_samples"]
        for pair in partition["patient_sample_map"]:
            observed_map[pair["patient_id"]] = pair["sample_id"]
    assert observed_map == patient_to_sample
    assert manifest["validation_results"]["passed"] is True


def test_manifest_records_schema_timestamp_sources_and_validation(
    split_project: dict[str, Path | list[dict[str, str]]],
) -> None:
    manifest = _generate(split_project)
    cohort_path = split_project["cohort"]
    matrix_manifest_path = split_project["matrix_manifest"]
    assert isinstance(cohort_path, Path) and isinstance(matrix_manifest_path, Path)

    assert manifest["schema_version"] == 1
    assert datetime.fromisoformat(manifest["created_at"]).utcoffset().total_seconds() == 8 * 3600
    assert manifest["random_seed"] == 20261006
    assert manifest["split_proportions"] == PROPORTIONS
    assert manifest["source_cohort"]["path"] == "data/processed/classification_cohort.tsv"
    assert (
        manifest["source_cohort"]["sha256"] == hashlib.sha256(cohort_path.read_bytes()).hexdigest()
    )
    assert manifest["expression_matrix_header_validation"]["matrix_manifest_sha256"] == (
        hashlib.sha256(matrix_manifest_path.read_bytes()).hexdigest()
    )
    assert manifest["validation_results"]["passed"] is True
    assert all(check["passed"] for check in manifest["validation_results"]["checks"])


def test_split_is_deterministic_for_fixed_seed_and_matrix_values_are_not_used(
    split_project: dict[str, Path | list[dict[str, str]]],
) -> None:
    first = _generate(split_project)
    matrix_path = split_project["matrix"]
    assert isinstance(matrix_path, Path)
    records = split_project["records"]
    assert isinstance(records, list)
    _write_matrix(matrix_path, [row["sample_id"] for row in records], values="-999999")
    second = _generate(split_project)

    assert first["partitions"] == second["partitions"]
    assert first["random_seed"] == second["random_seed"] == 20261006
    assert first["expression_matrix_header_validation"]["values_read_for_partitioning"] is False


def test_largest_remainder_rounding_is_explicit_and_reconciles() -> None:
    assert largest_remainder_sizes(844, PROPORTIONS, TIE_ORDER) == {
        "train": 506,
        "validation": 169,
        "test": 169,
    }
    assert largest_remainder_sizes(53, PROPORTIONS, TIE_ORDER) == {
        "train": 32,
        "validation": 11,
        "test": 10,
    }


def test_actual_sizes_and_pam50_counts_reconcile_with_cohort(
    split_project: dict[str, Path | list[dict[str, str]]],
) -> None:
    manifest = _generate(split_project)
    records = split_project["records"]
    assert isinstance(records, list)
    cohort_counts = Counter(row["pam50_normalized_label"] for row in records)
    assert manifest["cohort_pam50_counts"] == {
        label: cohort_counts[label] for label in PAM50_CLASS_ORDER
    }

    total_counts = Counter()
    for split, target in (("train", 30), ("validation", 10), ("test", 10)):
        partition = manifest["partitions"][split]
        assert partition["n_patients"] == target
        assert partition["n_samples"] == target
        total_counts.update(partition["pam50_counts"])
        assert all(partition["pam50_counts"][label] > 0 for label in PAM50_CLASS_ORDER)
        actual_fraction = target / len(records)
        assert abs(actual_fraction - PROPORTIONS[split]) <= 1 / len(records)
    assert total_counts == cohort_counts
    assert manifest["generation_method"]["stratified"] is True
    assert manifest["generation_method"]["expression_values_used_for_assignment"] is False


@pytest.mark.parametrize(
    ("label", "message"),
    [("", "Missing normalized PAM50 label"), ("Other subtype", "Invalid normalized PAM50 label")],
)
def test_missing_or_invalid_labels_fail_clearly(
    split_project: dict[str, Path | list[dict[str, str]]],
    label: str,
    message: str,
) -> None:
    records = split_project["records"]
    assert isinstance(records, list)
    changed = [dict(row) for row in records]
    changed[0]["pam50_normalized_label"] = label
    _write_cohort(split_project["cohort"], changed)  # type: ignore[arg-type]

    with pytest.raises(ValueError, match=message):
        _generate(split_project)


@pytest.mark.parametrize("duplicate_field", ["patient_id", "sample_id"])
def test_duplicate_identifiers_fail_clearly(
    split_project: dict[str, Path | list[dict[str, str]]],
    duplicate_field: str,
) -> None:
    records = split_project["records"]
    assert isinstance(records, list)
    changed = [dict(row) for row in records]
    changed[1][duplicate_field] = changed[0][duplicate_field]
    _write_cohort(split_project["cohort"], changed)  # type: ignore[arg-type]

    with pytest.raises(ValueError, match="duplicate"):
        _generate(split_project)


def test_expression_header_must_match_cohort_sample_ids(
    split_project: dict[str, Path | list[dict[str, str]]],
) -> None:
    matrix_path = split_project["matrix"]
    records = split_project["records"]
    assert isinstance(matrix_path, Path) and isinstance(records, list)
    _write_matrix(matrix_path, [row["sample_id"] for row in records[:-1]])

    with pytest.raises(ValueError, match="must match the approved cohort order exactly"):
        _generate(split_project)


def test_frozen_manifest_is_not_rewritten_or_reshuffled(
    split_project: dict[str, Path | list[dict[str, str]]],
) -> None:
    manifest = _generate(split_project)
    output_path = split_project["output"]
    assert isinstance(output_path, Path)
    assert write_frozen_split_manifest(manifest, output_path) is True
    frozen_bytes = output_path.read_bytes()

    same_assignment = _generate(split_project)
    assert write_frozen_split_manifest(same_assignment, output_path) is False
    assert output_path.read_bytes() == frozen_bytes

    changed_seed = _generate(split_project, seed=20261007)
    with pytest.raises(FileExistsError, match="Refusing to overwrite frozen split manifest"):
        write_frozen_split_manifest(changed_seed, output_path)
