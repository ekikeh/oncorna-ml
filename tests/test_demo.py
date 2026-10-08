from __future__ import annotations

import csv
import gzip
import hashlib
import json
from collections import Counter
from pathlib import Path
from typing import Any

import pytest

from oncorna.cohort import PAM50_CLASS_ORDER
from oncorna.demo import (
    build_demo_dataset,
    evenly_spaced_gene_indices,
    select_training_demo_samples,
    sha256_file,
)

COHORT_FIELDS = (
    "sample_id",
    "patient_id",
    "pam50_original_label",
    "pam50_normalized_label",
    "sample_type",
)
SEED = 20261006
EXPRESSION_SCALE = "log2(normalized_count + 1)"
ORIGINAL_LABELS = {
    "Luminal A": "LumA",
    "Luminal B": "LumB",
    "Basal-like": "Basal",
    "HER2-enriched": "Her2",
    "Normal-like": "Normal",
}


def _write_tsv(path: Path, fieldnames: tuple[str, ...], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def _write_matrix(
    path: Path, sample_ids: list[str], gene_ids: list[str]
) -> dict[str, dict[str, str]]:
    values: dict[str, dict[str, str]] = {}
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["gene_id", *sample_ids])
        for gene_index, gene_id in enumerate(gene_ids):
            row_values = {
                sample_id: f"{gene_index + sample_index / 10:.4f}"
                for sample_index, sample_id in enumerate(sample_ids)
            }
            values[gene_id] = row_values
            writer.writerow([gene_id, *[row_values[sample_id] for sample_id in sample_ids]])
    return values


@pytest.fixture
def demo_project(tmp_path: Path) -> dict[str, Any]:
    root = tmp_path
    records: list[dict[str, str]] = []
    partitions: dict[str, list[dict[str, str]]] = {
        "train": [],
        "validation": [],
        "test": [],
    }
    for class_index, normalized_label in enumerate(PAM50_CLASS_ORDER):
        for sample_index in range(10):
            index = class_index * 10 + sample_index
            row = {
                "sample_id": f"S{index:03d}",
                "patient_id": f"P{index:03d}",
                "pam50_original_label": ORIGINAL_LABELS[normalized_label],
                "pam50_normalized_label": normalized_label,
                "sample_type": "Primary Tumor",
            }
            records.append(row)
            if sample_index < 6:
                partition_name = "train"
            elif sample_index < 8:
                partition_name = "validation"
            else:
                partition_name = "test"
            partitions[partition_name].append(
                {"patient_id": row["patient_id"], "sample_id": row["sample_id"]}
            )

    cohort_path = root / "data/processed/classification_cohort.tsv"
    _write_tsv(cohort_path, COHORT_FIELDS, records)
    sample_ids = [row["sample_id"] for row in records]

    split_manifest_path = root / "data/processed/split_v1.json"
    split_payload: dict[str, Any] = {
        "schema_version": 1,
        "split_id": "split_v1",
        "split_unit": "patient_id",
        "random_seed": SEED,
        "validation_results": {"passed": True},
        "partitions": {},
    }
    for partition_name, mapping in partitions.items():
        patient_ids = [pair["patient_id"] for pair in mapping]
        partition_sample_ids = [pair["sample_id"] for pair in mapping]
        split_payload["partitions"][partition_name] = {
            "n_patients": len(mapping),
            "n_samples": len(mapping),
            "patient_ids": patient_ids,
            "sample_ids": partition_sample_ids,
            "patient_sample_map": mapping,
        }
    split_manifest_path.parent.mkdir(parents=True, exist_ok=True)
    split_manifest_path.write_text(
        json.dumps(split_payload, indent=2) + "\n",
        encoding="utf-8",
    )

    gene_ids = [f"G{index:03d}" for index in range(10)]
    gene_list_path = root / "data/processed/ml_gene_list.tsv"
    _write_tsv(gene_list_path, ("gene_id",), [{"gene_id": gene_id} for gene_id in gene_ids])

    expression_matrix_path = root / "data/processed/ml_expression_matrix.tsv.gz"
    source_values = _write_matrix(expression_matrix_path, sample_ids, gene_ids)
    matrix_manifest_path = root / "data/processed/ml_matrix_manifest.json"
    matrix_manifest = {
        "expression_units": EXPRESSION_SCALE,
        "sample_order": sample_ids,
        "gene_order": gene_ids,
        "matrix": {
            "path": "data/processed/ml_expression_matrix.tsv.gz",
            "sha256": sha256_file(expression_matrix_path),
            "shape": {
                "gene_rows": len(gene_ids),
                "sample_columns": len(sample_ids),
                "total_columns_including_gene_id": len(sample_ids) + 1,
            },
        },
        "gene_list": {
            "path": "data/processed/ml_gene_list.tsv",
            "sha256": sha256_file(gene_list_path),
        },
    }
    matrix_manifest_path.write_text(
        json.dumps(matrix_manifest, indent=2) + "\n",
        encoding="utf-8",
    )

    config_path = root / "configs/default.yaml"
    config_path.parent.mkdir(parents=True, exist_ok=True)
    config_path.write_text(f"random_seed: {SEED}\n", encoding="utf-8")
    split_payload["source_cohort"] = {
        "path": "data/processed/classification_cohort.tsv",
        "sha256": sha256_file(cohort_path),
        "patient_count": len(records),
        "sample_count": len(records),
    }
    split_payload["configuration"] = {
        "path": "configs/default.yaml",
        "sha256": sha256_file(config_path),
    }
    split_payload["expression_matrix_header_validation"] = {
        "path": "data/processed/ml_expression_matrix.tsv.gz",
        "matrix_manifest_path": "data/processed/ml_matrix_manifest.json",
        "matrix_manifest_sha256": sha256_file(matrix_manifest_path),
        "matrix_sha256_from_manifest": sha256_file(expression_matrix_path),
    }
    split_manifest_path.write_text(
        json.dumps(split_payload, indent=2) + "\n",
        encoding="utf-8",
    )
    demo_config_path = root / "configs/demo.yaml"
    demo_config_path.write_text(
        "partition: train\nsamples_per_class: 2\ngene_count: 4\n",
        encoding="utf-8",
    )

    return {
        "root": root,
        "cohort": cohort_path,
        "records": records,
        "split": split_manifest_path,
        "matrix": expression_matrix_path,
        "matrix_manifest": matrix_manifest_path,
        "gene_list": gene_list_path,
        "source_values": source_values,
        "gene_ids": gene_ids,
        "sample_ids": sample_ids,
        "config": config_path,
        "demo_config": demo_config_path,
        "expression_output": root / "data/demo/demo_expression.tsv",
        "metadata_output": root / "data/demo/demo_metadata.tsv",
        "manifest_output": root / "data/demo/demo_manifest.json",
    }


def _build(paths: dict[str, Any]):
    return build_demo_dataset(
        cohort_path=paths["cohort"],
        split_manifest_path=paths["split"],
        expression_matrix_path=paths["matrix"],
        gene_list_path=paths["gene_list"],
        matrix_manifest_path=paths["matrix_manifest"],
        project_config_path=paths["config"],
        demo_config_path=paths["demo_config"],
        project_root=paths["root"],
        expected_sample_count=50,
        expected_gene_count=10,
        expression_scale=EXPRESSION_SCALE,
        random_seed=SEED,
        samples_per_class=2,
        demo_gene_count=4,
        output_expression_path=paths["expression_output"],
        output_metadata_path=paths["metadata_output"],
        output_manifest_path=paths["manifest_output"],
    )


def _read_tsv(path: Path) -> tuple[list[str], list[list[str]]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        return next(reader), list(reader)


def test_evenly_spaced_gene_positions_are_fixed_and_ordered() -> None:
    assert evenly_spaced_gene_indices(10, 4) == [1, 3, 6, 8]
    assert evenly_spaced_gene_indices(17_623, 100) == sorted(
        evenly_spaced_gene_indices(17_623, 100)
    )
    assert len(set(evenly_spaced_gene_indices(17_623, 100))) == 100
    with pytest.raises(ValueError, match="Requested"):
        evenly_spaced_gene_indices(3, 4)


def test_demo_is_deterministic_and_has_expected_dimensions(demo_project: dict[str, Any]) -> None:
    paths = demo_project
    result = _build(paths)
    first_outputs = [
        paths["expression_output"].read_bytes(),
        paths["metadata_output"].read_bytes(),
        paths["manifest_output"].read_bytes(),
    ]
    second_result = _build(paths)
    second_outputs = [
        paths["expression_output"].read_bytes(),
        paths["metadata_output"].read_bytes(),
        paths["manifest_output"].read_bytes(),
    ]

    assert first_outputs == second_outputs
    assert result.sample_ids == second_result.sample_ids
    assert result.gene_ids == second_result.gene_ids
    assert len(result.sample_ids) == 10
    assert len(result.gene_ids) == 4

    expression_header, expression_rows = _read_tsv(paths["expression_output"])
    metadata_header, metadata_rows = _read_tsv(paths["metadata_output"])
    assert expression_header[0] == "gene_id"
    assert len(expression_rows) == 4
    assert len(expression_header) == 11
    assert metadata_header == [
        "sample_id",
        "patient_id",
        "pam50_original_label",
        "pam50_normalized_label",
        "split",
    ]
    assert len(metadata_rows) == 10
    assert len(expression_rows[0]) == 11


def test_metadata_alignment_train_partition_and_all_class_coverage(
    demo_project: dict[str, Any],
) -> None:
    paths = demo_project
    result = _build(paths)
    expression_header, _ = _read_tsv(paths["expression_output"])
    metadata_header, metadata_rows = _read_tsv(paths["metadata_output"])
    metadata = [dict(zip(metadata_header, row, strict=True)) for row in metadata_rows]
    cohort_by_sample = {row["sample_id"]: row for row in paths["records"]}

    assert expression_header[1:] == [row["sample_id"] for row in metadata]
    assert [row["sample_id"] for row in metadata] == list(result.sample_ids)
    assert all(row["split"] == "train" for row in metadata)
    assert all(row["sample_id"] in cohort_by_sample for row in metadata)
    assert all(
        row["patient_id"] == cohort_by_sample[row["sample_id"]]["patient_id"] for row in metadata
    )
    assert all(
        row["pam50_original_label"] == cohort_by_sample[row["sample_id"]]["pam50_original_label"]
        for row in metadata
    )
    assert all(
        row["pam50_normalized_label"]
        == cohort_by_sample[row["sample_id"]]["pam50_normalized_label"]
        for row in metadata
    )
    assert Counter(row["pam50_normalized_label"] for row in metadata) == {
        label: 2 for label in PAM50_CLASS_ORDER
    }
    assert result.class_counts == {label: 2 for label in PAM50_CLASS_ORDER}


def test_expression_tokens_are_copied_verbatim(demo_project: dict[str, Any]) -> None:
    paths = demo_project
    result = _build(paths)
    expression_header, expression_rows = _read_tsv(paths["expression_output"])
    selected_sample_ids = expression_header[1:]
    selected_gene_ids = [row[0] for row in expression_rows]

    assert selected_gene_ids == list(result.gene_ids)
    for row in expression_rows:
        gene_id = row[0]
        expected_tokens = [
            paths["source_values"][gene_id][sample_id] for sample_id in selected_sample_ids
        ]
        assert row[1:] == expected_tokens
    assert expression_rows


def test_selection_rules_ignore_expression_values_and_split_is_read_only(
    demo_project: dict[str, Any],
) -> None:
    paths = demo_project
    split_before = paths["split"].read_bytes()
    result = _build(paths)
    split_payload = json.loads(paths["split"].read_text(encoding="utf-8"))
    train_map = split_payload["partitions"]["train"]["patient_sample_map"]
    rows_with_arbitrary_values = [
        {**row, "expression_value_for_test_only": str(index * -7.25)}
        for index, row in enumerate(paths["records"])
    ]
    selected_with_arbitrary_values = select_training_demo_samples(
        rows_with_arbitrary_values,
        train_map,
        samples_per_class=2,
        random_seed=SEED,
    )
    selected_sample_ids = [row["sample_id"] for row in selected_with_arbitrary_values]
    selected_gene_ids = [paths["gene_ids"][index] for index in evenly_spaced_gene_indices(10, 4)]

    assert selected_sample_ids == list(result.sample_ids)
    assert selected_gene_ids == list(result.gene_ids)
    assert paths["split"].read_bytes() == split_before


def test_manifest_records_dimensions_rules_and_output_hashes(demo_project: dict[str, Any]) -> None:
    paths = demo_project
    result = _build(paths)
    manifest = json.loads(paths["manifest_output"].read_text(encoding="utf-8"))

    assert manifest["dimensions"] == {
        "gene_rows": 4,
        "sample_columns": 10,
        "total_expression_table_columns_including_gene_id": 11,
        "metadata_rows": 10,
    }
    assert manifest["sampling"]["partition"] == "train"
    assert manifest["sampling"]["pam50_counts"] == {label: 2 for label in PAM50_CLASS_ORDER}
    assert (
        manifest["gene_selection"]["position_formula"]
        == "floor(((2*j + 1) * N) / (2 * M)), j=0..M-1"
    )
    assert manifest["expression_values"]["additional_transformation"] == "none"
    assert (
        manifest["expression_values"]["copy_semantics"]
        == "selected source numeric tokens copied verbatim"
    )
    assert manifest["outputs"]["expression"]["sha256"] == result.expression_sha256
    assert manifest["outputs"]["metadata"]["sha256"] == result.metadata_sha256
    assert (
        result.manifest_sha256 == hashlib.sha256(paths["manifest_output"].read_bytes()).hexdigest()
    )
