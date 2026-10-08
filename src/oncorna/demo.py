"""Build a compact, reproducible Phase 1G demo from frozen Phase 1C/1E/1F data.

Samples are selected from the frozen training partition by a SHA-256 rank over
patient IDs and the configured seed. Genes are selected by evenly spaced row
positions in the frozen Phase 1E gene list. Expression values are copied
verbatim; no value-based selection or transformation occurs here.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence

from oncorna.cohort import PAM50_CLASS_ORDER
from oncorna.qc import read_classification_cohort

DEMO_METADATA_FIELDS = (
    "sample_id",
    "patient_id",
    "pam50_original_label",
    "pam50_normalized_label",
    "split",
)
MISSING_LABELS = {"", "--", "n/a", "na", "nan", "none", "not available", "null"}
DEMO_HASH_NAMESPACE = "oncorna-ml-phase-1g-demo-v1"


@dataclass(frozen=True)
class DemoBuildResult:
    """Summary of the committed demo artifacts generated locally."""

    sample_ids: tuple[str, ...]
    patient_ids: tuple[str, ...]
    gene_ids: tuple[str, ...]
    class_counts: dict[str, int]
    expression_path: Path
    metadata_path: Path
    manifest_path: Path
    expression_sha256: str
    metadata_sha256: str
    manifest_sha256: str


def sha256_file(path: Path, *, chunk_size: int = 1024 * 1024) -> str:
    """Return a streaming SHA-256 digest for a file."""
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _project_relative(path: Path, project_root: Path) -> str:
    try:
        return path.relative_to(project_root).as_posix()
    except ValueError as exc:
        raise ValueError(f"Demo input/output must be inside the project: {path}") from exc


def _input_record(path: Path, project_root: Path) -> dict[str, str]:
    return {
        "path": _project_relative(path, project_root),
        "sha256": sha256_file(path),
    }


def evenly_spaced_gene_indices(total_genes: int, demo_gene_count: int) -> list[int]:
    """Return deterministic, evenly spaced zero-based row positions.

    Position ``j`` is ``floor(((2*j + 1) * N) / (2*M))`` where ``N`` is the
    source-list length and ``M`` the requested demo-gene count. This samples
    across the list order without using expression values or labels.
    """
    if total_genes < 1 or demo_gene_count < 1:
        raise ValueError("Total and demo gene counts must be positive")
    if demo_gene_count > total_genes:
        raise ValueError(
            f"Requested {demo_gene_count} demo genes from a {total_genes}-gene source list"
        )
    indices = [
        ((2 * index + 1) * total_genes) // (2 * demo_gene_count) for index in range(demo_gene_count)
    ]
    if len(indices) != len(set(indices)):
        raise ValueError("Evenly spaced gene-position rule produced duplicate indices")
    return indices


def _read_gene_list(path: Path, *, expected_gene_count: int) -> list[str]:
    with Path(path).open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if "gene_id" not in (reader.fieldnames or ()):
            raise ValueError(f"Gene list {path} must contain a gene_id column")
        gene_ids = [(row.get("gene_id") or "").strip() for row in reader]
    if len(gene_ids) != expected_gene_count:
        raise ValueError(f"Gene list has {len(gene_ids)} rows; expected {expected_gene_count}")
    if any(not gene_id for gene_id in gene_ids):
        raise ValueError("Gene list contains a blank gene_id")
    if len(gene_ids) != len(set(gene_ids)):
        raise ValueError("Gene list contains duplicate gene IDs")
    return gene_ids


def _validate_cohort_labels(
    cohort_path: Path,
    *,
    expected_sample_count: int,
) -> list[dict[str, str]]:
    records = read_classification_cohort(
        cohort_path,
        expected_sample_count=expected_sample_count,
    )
    by_patient: dict[str, dict[str, str]] = {}
    sample_ids: set[str] = set()
    for row in records:
        patient_id = row["patient_id"]
        sample_id = row["sample_id"]
        if patient_id in by_patient:
            raise ValueError(f"Duplicate cohort patient_id: {patient_id}")
        if sample_id in sample_ids:
            raise ValueError(f"Duplicate cohort sample_id: {sample_id}")
        label = row.get("pam50_normalized_label", "").strip()
        if label.casefold() in MISSING_LABELS:
            raise ValueError(f"Missing normalized PAM50 label for sample {sample_id}")
        if label not in PAM50_CLASS_ORDER:
            raise ValueError(f"Invalid normalized PAM50 label {label!r} for sample {sample_id}")
        by_patient[patient_id] = row
        sample_ids.add(sample_id)
    return records


def _validate_frozen_split(
    split_manifest: Mapping[str, Any],
    cohort_records: Sequence[Mapping[str, str]],
    *,
    random_seed: int,
) -> list[dict[str, str]]:
    """Validate the immutable split manifest and return its training patient/sample map."""
    if split_manifest.get("schema_version") != 1 or split_manifest.get("split_id") != "split_v1":
        raise ValueError("Demo generation requires the frozen split_v1 manifest")
    if split_manifest.get("split_unit") != "patient_id":
        raise ValueError("Frozen split_v1 must use patient_id as its split unit")
    if split_manifest.get("random_seed") != random_seed:
        raise ValueError("Demo selection seed must match the seed recorded in frozen split_v1")
    if split_manifest.get("validation_results", {}).get("passed") is not True:
        raise ValueError("Frozen split_v1 does not contain passing validation results")

    cohort_map = {(row["patient_id"], row["sample_id"]) for row in cohort_records}
    observed_pairs: list[tuple[str, str]] = []
    train_pairs: list[dict[str, str]] = []
    for partition_name in ("train", "validation", "test"):
        partition = split_manifest.get("partitions", {}).get(partition_name)
        if not isinstance(partition, dict):
            raise ValueError(f"Frozen split_v1 is missing partition {partition_name!r}")
        patient_ids = partition.get("patient_ids", [])
        sample_ids = partition.get("sample_ids", [])
        mapping = partition.get("patient_sample_map", [])
        if len(patient_ids) != len(set(patient_ids)) or len(sample_ids) != len(set(sample_ids)):
            raise ValueError(f"Frozen split_v1 has duplicate IDs in {partition_name}")
        mapped_pairs = [(pair.get("patient_id", ""), pair.get("sample_id", "")) for pair in mapping]
        if len(mapped_pairs) != len(mapping) or len(mapped_pairs) != len(patient_ids):
            raise ValueError(f"Frozen split_v1 mapping is incomplete in {partition_name}")
        if [pair[0] for pair in mapped_pairs] != patient_ids:
            raise ValueError(
                f"Frozen split_v1 patient ordering/mapping differs in {partition_name}"
            )
        if [pair[1] for pair in mapped_pairs] != sample_ids:
            raise ValueError(f"Frozen split_v1 sample ordering/mapping differs in {partition_name}")
        if partition.get("n_patients") != len(patient_ids) or partition.get("n_samples") != len(
            sample_ids
        ):
            raise ValueError(f"Frozen split_v1 counts disagree in {partition_name}")
        if any(pair not in cohort_map for pair in mapped_pairs):
            raise ValueError(
                f"Frozen split_v1 contains an unapproved patient/sample pair in {partition_name}"
            )
        observed_pairs.extend(mapped_pairs)
        if partition_name == "train":
            train_pairs = [
                {"patient_id": patient_id, "sample_id": sample_id}
                for patient_id, sample_id in mapped_pairs
            ]

    if len(observed_pairs) != len(set(observed_pairs)):
        raise ValueError("Frozen split_v1 repeats a patient/sample pair across partitions")
    if set(observed_pairs) != cohort_map:
        raise ValueError("Frozen split_v1 does not partition the approved cohort exactly")
    return train_pairs


def select_training_demo_samples(
    cohort_records: Sequence[Mapping[str, str]],
    train_patient_sample_map: Sequence[Mapping[str, str]],
    *,
    samples_per_class: int,
    random_seed: int,
) -> list[dict[str, str]]:
    """Select a balanced demo subset from training patients using only IDs and labels.

    Within each normalized PAM50 class, rank training patients by SHA-256 of
    ``oncorna-ml-phase-1g-demo-v1|<seed>|<patient_id>`` and take the first K.
    This deterministic hash ranking is independent of expression values.
    """
    if samples_per_class < 1:
        raise ValueError("samples_per_class must be positive")
    by_patient = {row["patient_id"]: row for row in cohort_records}
    train_rows: list[dict[str, str]] = []
    seen_patients: set[str] = set()
    seen_samples: set[str] = set()
    for pair in train_patient_sample_map:
        patient_id = pair.get("patient_id", "")
        sample_id = pair.get("sample_id", "")
        if patient_id in seen_patients or sample_id in seen_samples:
            raise ValueError("Training patient/sample mapping contains duplicate identifiers")
        if patient_id not in by_patient or by_patient[patient_id]["sample_id"] != sample_id:
            raise ValueError(
                f"Training mapping is not present in approved cohort: {patient_id}/{sample_id}"
            )
        seen_patients.add(patient_id)
        seen_samples.add(sample_id)
        train_rows.append(by_patient[patient_id])

    selected: list[dict[str, str]] = []
    for label in PAM50_CLASS_ORDER:
        candidates = [row for row in train_rows if row["pam50_normalized_label"].strip() == label]
        if len(candidates) < samples_per_class:
            raise ValueError(
                f"Training partition has only {len(candidates)} {label} samples; "
                f"cannot select {samples_per_class}"
            )
        ranked = sorted(
            candidates,
            key=lambda row: (
                hashlib.sha256(
                    f"{DEMO_HASH_NAMESPACE}|{random_seed}|{row['patient_id']}".encode("utf-8")
                ).hexdigest(),
                row["patient_id"],
            ),
        )
        selected.extend(dict(row) for row in ranked[:samples_per_class])
    return selected


def _write_tsv_temp(
    path: Path,
    rows: Sequence[Mapping[str, Any]],
    fields: Sequence[str],
) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(fields),
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


def _write_expression_temp(
    *,
    source_matrix_path: Path,
    temp_output_path: Path,
    cohort_sample_ids: Sequence[str],
    demo_sample_ids: Sequence[str],
    selected_gene_ids: Sequence[str],
    expected_gene_count: int,
) -> None:
    """Copy exactly the selected source tokens into a small uncompressed TSV."""
    with gzip.open(source_matrix_path, "rt", encoding="utf-8-sig", newline="") as source:
        reader = csv.reader(source, delimiter="\t")
        header = next(reader, None)
        expected_header = ["gene_id", *cohort_sample_ids]
        if header != expected_header:
            raise ValueError("Full ML matrix header differs from the approved cohort order")
        column_by_sample = {
            sample_id: index + 1 for index, sample_id in enumerate(cohort_sample_ids)
        }
        if any(sample_id not in column_by_sample for sample_id in demo_sample_ids):
            raise ValueError("Selected demo sample is absent from the full expression matrix")
        sample_indices = [column_by_sample[sample_id] for sample_id in demo_sample_ids]
        selected_gene_set = set(selected_gene_ids)
        next_gene_index = 0
        seen_gene_ids: set[str] = set()
        source_gene_count = 0
        with temp_output_path.open("w", encoding="utf-8", newline="") as output:
            writer = csv.writer(output, delimiter="\t", lineterminator="\n")
            writer.writerow(["gene_id", *demo_sample_ids])
            for row_number, row in enumerate(reader, start=2):
                if len(row) != len(expected_header):
                    raise ValueError(
                        f"Full ML matrix row {row_number} has {len(row)} fields; "
                        f"expected {len(expected_header)}"
                    )
                gene_id = row[0]
                if not gene_id or gene_id in seen_gene_ids:
                    raise ValueError(
                        f"Full ML matrix has a blank or duplicate gene ID: {gene_id!r}"
                    )
                seen_gene_ids.add(gene_id)
                source_gene_count += 1
                if gene_id not in selected_gene_set:
                    continue
                if (
                    next_gene_index >= len(selected_gene_ids)
                    or gene_id != selected_gene_ids[next_gene_index]
                ):
                    raise ValueError(
                        "Selected gene order differs from the full matrix source order"
                    )
                values: list[str] = []
                for sample_id, column_index in zip(demo_sample_ids, sample_indices, strict=True):
                    source_token = row[column_index]
                    try:
                        value = float(source_token.strip())
                    except ValueError as exc:
                        raise ValueError(
                            f"Non-numeric source value for {gene_id}/{sample_id}: {source_token!r}"
                        ) from exc
                    if not math.isfinite(value):
                        raise ValueError(
                            f"Non-finite source value for {gene_id}/{sample_id}: {source_token!r}"
                        )
                    values.append(source_token)
                writer.writerow([gene_id, *values])
                next_gene_index += 1

    if source_gene_count != expected_gene_count:
        raise ValueError(
            f"Full ML matrix has {source_gene_count} genes; expected {expected_gene_count}"
        )
    if next_gene_index != len(selected_gene_ids):
        missing = selected_gene_ids[next_gene_index : next_gene_index + 5]
        raise ValueError(f"Selected demo genes are missing from matrix: {missing}")


def build_demo_dataset(
    *,
    cohort_path: Path,
    split_manifest_path: Path,
    expression_matrix_path: Path,
    gene_list_path: Path,
    matrix_manifest_path: Path,
    project_config_path: Path,
    demo_config_path: Path,
    project_root: Path,
    expected_sample_count: int,
    expected_gene_count: int,
    expression_scale: str,
    random_seed: int,
    samples_per_class: int,
    demo_gene_count: int,
    output_expression_path: Path,
    output_metadata_path: Path,
    output_manifest_path: Path,
) -> DemoBuildResult:
    """Create deterministic demo tables from existing Phase 1C/1E/1F outputs."""
    project_root = Path(project_root).resolve()
    cohort_path = Path(cohort_path).resolve()
    split_manifest_path = Path(split_manifest_path).resolve()
    expression_matrix_path = Path(expression_matrix_path).resolve()
    gene_list_path = Path(gene_list_path).resolve()
    matrix_manifest_path = Path(matrix_manifest_path).resolve()
    project_config_path = Path(project_config_path).resolve()
    demo_config_path = Path(demo_config_path).resolve()
    output_expression_path = Path(output_expression_path).resolve()
    output_metadata_path = Path(output_metadata_path).resolve()
    output_manifest_path = Path(output_manifest_path).resolve()

    demo_dir = (project_root / "data" / "demo").resolve()
    output_paths = (output_expression_path, output_metadata_path, output_manifest_path)
    if any(not path.is_relative_to(demo_dir) for path in output_paths):
        raise ValueError("Phase 1G outputs must remain under data/demo/")
    if len(set(output_paths)) != len(output_paths):
        raise ValueError("Demo output paths must be distinct")
    if expression_scale.strip() != "log2(normalized_count + 1)":
        raise ValueError("Phase 1G requires expression units log2(normalized_count + 1)")
    if expected_sample_count < 1 or expected_gene_count < 1 or random_seed < 0:
        raise ValueError("Expected dimensions and random seed must be nonnegative/positive")

    split_sha_before = sha256_file(split_manifest_path)
    split_manifest = json.loads(split_manifest_path.read_text(encoding="utf-8"))
    cohort_records = _validate_cohort_labels(
        cohort_path,
        expected_sample_count=expected_sample_count,
    )
    if split_manifest.get("source_cohort", {}).get("path") != _project_relative(
        cohort_path, project_root
    ):
        raise ValueError("Frozen split_v1 was generated from a different cohort path")
    if split_manifest.get("source_cohort", {}).get("sha256") != sha256_file(cohort_path):
        raise ValueError("Frozen split_v1 cohort checksum does not match the approved cohort")
    if split_manifest.get("configuration", {}).get("path") != _project_relative(
        project_config_path, project_root
    ):
        raise ValueError("Frozen split_v1 was generated from a different project configuration")
    if split_manifest.get("configuration", {}).get("sha256") != sha256_file(project_config_path):
        raise ValueError("Frozen split_v1 project configuration checksum does not match")
    split_matrix_check = split_manifest.get("expression_matrix_header_validation", {})
    if split_matrix_check.get("path") != _project_relative(expression_matrix_path, project_root):
        raise ValueError("Frozen split_v1 validated a different expression matrix")
    if split_matrix_check.get("matrix_manifest_path") != _project_relative(
        matrix_manifest_path, project_root
    ):
        raise ValueError("Frozen split_v1 validated a different Phase 1E matrix manifest")
    if split_matrix_check.get("matrix_manifest_sha256") != sha256_file(matrix_manifest_path):
        raise ValueError("Frozen split_v1 matrix-manifest checksum does not match")
    cohort_sample_ids = [row["sample_id"] for row in cohort_records]
    cohort_patient_pairs = {(row["patient_id"], row["sample_id"]) for row in cohort_records}
    train_map = _validate_frozen_split(
        split_manifest,
        cohort_records,
        random_seed=random_seed,
    )
    selected_sample_records = select_training_demo_samples(
        cohort_records,
        train_map,
        samples_per_class=samples_per_class,
        random_seed=random_seed,
    )
    selected_sample_ids = [row["sample_id"] for row in selected_sample_records]
    selected_patient_ids = [row["patient_id"] for row in selected_sample_records]
    if len(selected_sample_ids) != len(set(selected_sample_ids)):
        raise ValueError("Demo sampling selected duplicate sample IDs")
    if len(selected_patient_ids) != len(set(selected_patient_ids)):
        raise ValueError("Demo sampling selected duplicate patient IDs")
    if not set(zip(selected_patient_ids, selected_sample_ids, strict=True)) <= cohort_patient_pairs:
        raise ValueError(
            "Demo sampling selected a patient/sample mapping outside the approved cohort"
        )

    matrix_manifest = json.loads(matrix_manifest_path.read_text(encoding="utf-8"))
    if matrix_manifest.get("expression_units") != expression_scale:
        raise ValueError("Phase 1E manifest expression units do not match configuration")
    matrix_shape = matrix_manifest.get("matrix", {}).get("shape", {})
    if matrix_shape.get("sample_columns") != expected_sample_count:
        raise ValueError("Phase 1E manifest sample count does not match configuration")
    if matrix_shape.get("gene_rows") != expected_gene_count:
        raise ValueError("Phase 1E manifest gene count does not match configuration")
    if matrix_manifest.get("sample_order") != cohort_sample_ids:
        raise ValueError("Phase 1E manifest sample order differs from classification cohort")
    if matrix_manifest.get("matrix", {}).get("path") != _project_relative(
        expression_matrix_path, project_root
    ):
        raise ValueError("Phase 1E manifest points to a different expression matrix")
    if matrix_manifest.get("gene_list", {}).get("path") != _project_relative(
        gene_list_path, project_root
    ):
        raise ValueError("Phase 1E manifest points to a different retained-gene list")

    gene_ids = _read_gene_list(gene_list_path, expected_gene_count=expected_gene_count)
    if matrix_manifest.get("gene_order") != gene_ids:
        raise ValueError("Phase 1E manifest gene order differs from the retained-gene list")
    selected_gene_indices = evenly_spaced_gene_indices(len(gene_ids), demo_gene_count)
    selected_gene_ids = [gene_ids[index] for index in selected_gene_indices]
    if len(selected_gene_ids) != len(set(selected_gene_ids)):
        raise ValueError("Demo gene selection contains duplicate gene IDs")

    matrix_digest = sha256_file(expression_matrix_path)
    gene_list_digest = sha256_file(gene_list_path)
    if matrix_manifest.get("matrix", {}).get("sha256") != matrix_digest:
        raise ValueError("Full expression matrix does not match Phase 1E manifest checksum")
    if split_matrix_check.get("matrix_sha256_from_manifest") != matrix_digest:
        raise ValueError("Frozen split_v1 records a different Phase 1E matrix checksum")
    if matrix_manifest.get("gene_list", {}).get("sha256") != gene_list_digest:
        raise ValueError("Retained-gene list does not match Phase 1E manifest checksum")

    output_expression_path.parent.mkdir(parents=True, exist_ok=True)
    expression_temp = output_expression_path.with_name(output_expression_path.name + ".tmp")
    metadata_temp = output_metadata_path.with_name(output_metadata_path.name + ".tmp")
    manifest_temp = output_manifest_path.with_name(output_manifest_path.name + ".tmp")
    temp_paths = (expression_temp, metadata_temp, manifest_temp)
    try:
        _write_expression_temp(
            source_matrix_path=expression_matrix_path,
            temp_output_path=expression_temp,
            cohort_sample_ids=cohort_sample_ids,
            demo_sample_ids=selected_sample_ids,
            selected_gene_ids=selected_gene_ids,
            expected_gene_count=expected_gene_count,
        )
        metadata_rows = [
            {
                "sample_id": row["sample_id"],
                "patient_id": row["patient_id"],
                "pam50_original_label": row["pam50_original_label"],
                "pam50_normalized_label": row["pam50_normalized_label"],
                "split": "train",
            }
            for row in selected_sample_records
        ]
        _write_tsv_temp(metadata_temp, metadata_rows, DEMO_METADATA_FIELDS)
        class_counts = {
            label: sum(row["pam50_normalized_label"] == label for row in metadata_rows)
            for label in PAM50_CLASS_ORDER
        }
        if any(class_counts[label] != samples_per_class for label in PAM50_CLASS_ORDER):
            raise ValueError("Demo sample selection did not cover every configured PAM50 class")

        input_files = {
            "classification_cohort": _input_record(cohort_path, project_root),
            "frozen_split_v1": _input_record(split_manifest_path, project_root),
            "full_expression_matrix": {
                **_input_record(expression_matrix_path, project_root),
                "sha256_declared_by_matrix_manifest": matrix_manifest["matrix"].get("sha256"),
            },
            "retained_gene_list": {
                **_input_record(gene_list_path, project_root),
                "sha256_declared_by_matrix_manifest": matrix_manifest["gene_list"].get("sha256"),
            },
            "matrix_manifest": _input_record(matrix_manifest_path, project_root),
            "project_configuration": _input_record(project_config_path, project_root),
            "demo_configuration": _input_record(demo_config_path, project_root),
        }
        expression_sha = sha256_file(expression_temp)
        metadata_sha = sha256_file(metadata_temp)
        demo_manifest = {
            "schema_version": 1,
            "dataset_id": "oncorna_ml_phase1g_demo_v1",
            "purpose": (
                "Small, class-balanced training-partition sample for inspecting table structure; "
                "not for biological feature selection or performance estimation."
            ),
            "inputs": input_files,
            "sampling": {
                "partition": "train",
                "method": (
                    "Within each normalized PAM50 class, rank training patients by SHA-256 of "
                    f"{DEMO_HASH_NAMESPACE}|{random_seed}|<patient_id>; choose the first "
                    f"{samples_per_class} in rank order."
                ),
                "hash_namespace": DEMO_HASH_NAMESPACE,
                "seed": random_seed,
                "pam50_class_order": list(PAM50_CLASS_ORDER),
                "samples_per_class": samples_per_class,
                "sample_count": len(selected_sample_ids),
                "pam50_counts": class_counts,
                "selected_sample_ids": selected_sample_ids,
            },
            "gene_selection": {
                "source_gene_list": _project_relative(gene_list_path, project_root),
                "source_gene_count": len(gene_ids),
                "demo_gene_count": len(selected_gene_ids),
                "method": "evenly spaced zero-based positions in Phase 1E gene-list order",
                "position_formula": "floor(((2*j + 1) * N) / (2 * M)), j=0..M-1",
                "selected_zero_based_positions": selected_gene_indices,
                "selected_gene_ids": selected_gene_ids,
                "biological_feature_selection": False,
            },
            "expression_values": {
                "units": expression_scale,
                "additional_transformation": "none",
                "copy_semantics": "selected source numeric tokens copied verbatim",
            },
            "dimensions": {
                "gene_rows": len(selected_gene_ids),
                "sample_columns": len(selected_sample_ids),
                "total_expression_table_columns_including_gene_id": len(selected_sample_ids) + 1,
                "metadata_rows": len(metadata_rows),
            },
            "outputs": {
                "expression": {
                    "path": _project_relative(output_expression_path, project_root),
                    "format": "UTF-8 TSV",
                    "sha256": expression_sha,
                },
                "metadata": {
                    "path": _project_relative(output_metadata_path, project_root),
                    "format": "UTF-8 TSV",
                    "sha256": metadata_sha,
                },
            },
        }
        with manifest_temp.open("w", encoding="utf-8", newline="\n") as handle:
            json.dump(demo_manifest, handle, indent=2, sort_keys=True)
            handle.write("\n")

        if sha256_file(split_manifest_path) != split_sha_before:
            raise RuntimeError("Frozen split manifest changed during demo generation")
        expression_temp.replace(output_expression_path)
        metadata_temp.replace(output_metadata_path)
        manifest_temp.replace(output_manifest_path)
    finally:
        for path in temp_paths:
            if path.exists():
                path.unlink()

    return DemoBuildResult(
        sample_ids=tuple(selected_sample_ids),
        patient_ids=tuple(selected_patient_ids),
        gene_ids=tuple(selected_gene_ids),
        class_counts=class_counts,
        expression_path=output_expression_path,
        metadata_path=output_metadata_path,
        manifest_path=output_manifest_path,
        expression_sha256=sha256_file(output_expression_path),
        metadata_sha256=sha256_file(output_metadata_path),
        manifest_sha256=sha256_file(output_manifest_path),
    )
