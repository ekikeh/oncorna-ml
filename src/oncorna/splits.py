"""Deterministic, patient-level train/validation/test partitioning.

Only cohort identifiers and normalized PAM50 labels drive assignment. The
expression matrix is consulted header-only to verify sample membership; its
expression values are never loaded or used for splitting.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
from collections import Counter
from datetime import datetime
from pathlib import Path
from typing import Any, Mapping, Sequence
from zoneinfo import ZoneInfo

import sklearn
from sklearn.model_selection import train_test_split

from oncorna.cohort import PAM50_CLASS_ORDER
from oncorna.qc import read_classification_cohort

PARTITION_ORDER = ("train", "validation", "test")
MISSING_LABELS = {"", "--", "n/a", "na", "nan", "none", "not available", "null"}


class StratificationError(ValueError):
    """Raised when the requested split sizes cannot support label stratification."""


def sha256_file(path: Path) -> str:
    """Return a streaming SHA-256 digest for a file."""
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _manifest_path(path: Path, project_root: Path) -> str:
    try:
        return path.relative_to(project_root).as_posix()
    except ValueError:
        return path.as_posix()


def largest_remainder_sizes(
    total: int,
    proportions: Mapping[str, float],
    tie_break_order: Sequence[str],
) -> dict[str, int]:
    """Apportion integer partition sizes with Hamilton's largest-remainder rule.

    Floors are assigned first. Remaining patients go to the largest fractional
    remainders; ties follow ``tie_break_order``. For 844 patients at 60/20/20,
    this yields 506/169/169.
    """
    keys = tuple(proportions)
    if total < 1:
        raise ValueError("Total patient count must be positive")
    if set(keys) != set(PARTITION_ORDER) or len(keys) != len(PARTITION_ORDER):
        raise ValueError(f"Proportions must define exactly {PARTITION_ORDER}")
    if len(tie_break_order) != len(keys) or set(tie_break_order) != set(keys):
        raise ValueError("Rounding tie-break order must list every partition exactly once")
    if any(not math.isfinite(float(proportions[key])) or proportions[key] <= 0 for key in keys):
        raise ValueError("Every split proportion must be finite and greater than zero")
    if not math.isclose(sum(proportions.values()), 1.0, rel_tol=0, abs_tol=1e-12):
        raise ValueError("Split proportions must sum to 1.0")

    exact_sizes = {key: total * float(proportions[key]) for key in keys}
    sizes = {key: math.floor(value) for key, value in exact_sizes.items()}
    remaining = total - sum(sizes.values())
    tie_rank = {key: index for index, key in enumerate(tie_break_order)}
    ranked = sorted(
        keys,
        key=lambda key: (-(exact_sizes[key] - sizes[key]), tie_rank[key]),
    )
    for key in ranked[:remaining]:
        sizes[key] += 1
    if sum(sizes.values()) != total:
        raise RuntimeError("Largest-remainder apportionment failed to preserve the total")
    return {key: sizes[key] for key in PARTITION_ORDER}


def _read_validated_cohort(
    cohort_path: Path,
    *,
    expected_patient_count: int,
) -> tuple[list[dict[str, str]], list[str]]:
    records = read_classification_cohort(
        cohort_path,
        expected_sample_count=expected_patient_count,
    )
    patient_ids = [row["patient_id"] for row in records]
    sample_ids = [row["sample_id"] for row in records]
    if len(patient_ids) != len(set(patient_ids)):
        raise ValueError("Classification cohort contains duplicate patient IDs")
    if len(sample_ids) != len(set(sample_ids)):
        raise ValueError("Classification cohort contains duplicate sample IDs")
    if len(patient_ids) != len(sample_ids):
        raise ValueError(
            "Patient-level split requires exactly one approved expression sample per patient"
        )

    labels: list[str] = []
    for row in records:
        sample_id = row["sample_id"]
        label = row.get("pam50_normalized_label", "").strip()
        if label.casefold() in MISSING_LABELS:
            raise ValueError(f"Missing normalized PAM50 label for sample {sample_id}")
        if label not in PAM50_CLASS_ORDER:
            raise ValueError(
                f"Invalid normalized PAM50 label {label!r} for sample {sample_id}; "
                f"expected one of {list(PAM50_CLASS_ORDER)}"
            )
        labels.append(label)
    observed_classes = set(labels)
    missing_classes = [label for label in PAM50_CLASS_ORDER if label not in observed_classes]
    if missing_classes:
        raise ValueError(f"Classification cohort is missing PAM50 classes: {missing_classes}")
    return records, labels


def _read_expression_matrix_header(
    expression_matrix_path: Path,
    matrix_manifest_path: Path,
    *,
    project_root: Path,
    cohort_sample_ids: Sequence[str],
    expected_gene_count: int,
    expression_scale: str,
) -> dict[str, Any]:
    """Validate matrix sample IDs using its header only; never read expression rows."""
    with gzip.open(expression_matrix_path, "rt", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, None)
    if not header or len(header) < 2 or header[0] != "gene_id":
        raise ValueError("ML expression matrix must start with a gene_id column and sample IDs")
    matrix_sample_ids = header[1:]
    if any(not sample_id for sample_id in matrix_sample_ids):
        raise ValueError("ML expression matrix header contains a blank sample ID")
    if len(matrix_sample_ids) != len(set(matrix_sample_ids)):
        raise ValueError("ML expression matrix header contains duplicate sample IDs")
    if list(matrix_sample_ids) != list(cohort_sample_ids):
        missing = sorted(set(cohort_sample_ids) - set(matrix_sample_ids))
        extra = sorted(set(matrix_sample_ids) - set(cohort_sample_ids))
        raise ValueError(
            "ML expression matrix sample header must match the approved cohort order exactly: "
            f"missing={missing[:5]}, extra={extra[:5]}"
        )

    with Path(matrix_manifest_path).open("r", encoding="utf-8") as handle:
        matrix_manifest = json.load(handle)
    expected_matrix_path = _manifest_path(Path(expression_matrix_path).resolve(), project_root)
    declared_matrix_path = matrix_manifest.get("matrix", {}).get("path")
    if declared_matrix_path != expected_matrix_path:
        raise ValueError(
            "ML matrix manifest path does not identify the configured expression matrix: "
            f"{declared_matrix_path!r} != {expected_matrix_path!r}"
        )
    manifest_sample_ids = matrix_manifest.get("sample_order")
    if manifest_sample_ids != list(cohort_sample_ids):
        raise ValueError("ML matrix manifest sample order does not match the approved cohort")
    shape = matrix_manifest.get("matrix", {}).get("shape", {})
    if shape.get("sample_columns") != len(cohort_sample_ids):
        raise ValueError("ML matrix manifest sample dimension does not match the cohort")
    if shape.get("total_columns_including_gene_id") != len(cohort_sample_ids) + 1:
        raise ValueError("ML matrix manifest total column count does not match its sample header")
    if shape.get("gene_rows") != expected_gene_count:
        raise ValueError(
            f"ML matrix manifest records {shape.get('gene_rows')} genes; "
            f"expected {expected_gene_count}"
        )
    if matrix_manifest.get("expression_units") != expression_scale:
        raise ValueError("ML matrix manifest expression units do not match project configuration")

    return {
        "sample_column_count": len(matrix_sample_ids),
        "gene_row_count": shape["gene_rows"],
        "sample_order_matches_cohort": True,
        "values_read_for_partitioning": False,
        "matrix_sha256_from_manifest": matrix_manifest.get("matrix", {}).get("sha256"),
    }


def _stratified_partition_indices(
    labels: Sequence[str],
    *,
    target_sizes: Mapping[str, int],
    random_seed: int,
    test_seed_offset: int,
    validation_seed_offset: int,
) -> dict[str, list[int]]:
    """Split indices in two deterministic stratified stages."""
    if min(test_seed_offset, validation_seed_offset) < 0:
        raise ValueError("Random-state offsets must be nonnegative")
    test_random_state = random_seed + test_seed_offset
    validation_random_state = random_seed + validation_seed_offset
    max_random_state = 2**32 - 1
    if max(test_random_state, validation_random_state) > max_random_state:
        raise ValueError("Random seed plus offsets must be less than 2**32")

    indices = list(range(len(labels)))
    try:
        train_validation_indices, test_indices = train_test_split(
            indices,
            test_size=target_sizes["test"],
            random_state=test_random_state,
            shuffle=True,
            stratify=list(labels),
        )
        train_indices, validation_indices = train_test_split(
            train_validation_indices,
            test_size=target_sizes["validation"],
            random_state=validation_random_state,
            shuffle=True,
            stratify=[labels[index] for index in train_validation_indices],
        )
    except ValueError as exc:
        counts = dict(Counter(labels))
        raise StratificationError(
            "Stratified train/validation/test split is infeasible for the requested sizes "
            f"and PAM50 counts {counts}: {exc}"
        ) from exc

    # Assignment is complete; serialize members in frozen cohort row order.
    return {
        "train": sorted(int(index) for index in train_indices),
        "validation": sorted(int(index) for index in validation_indices),
        "test": sorted(int(index) for index in test_indices),
    }


def _partition_records(
    records: Sequence[Mapping[str, str]],
    labels: Sequence[str],
    indices_by_partition: Mapping[str, Sequence[int]],
) -> dict[str, dict[str, Any]]:
    partitions: dict[str, dict[str, Any]] = {}
    for partition in PARTITION_ORDER:
        selected_records = [records[index] for index in indices_by_partition[partition]]
        counts = Counter(labels[index] for index in indices_by_partition[partition])
        patient_ids = [row["patient_id"] for row in selected_records]
        sample_ids = [row["sample_id"] for row in selected_records]
        partitions[partition] = {
            "n_patients": len(selected_records),
            "n_samples": len(selected_records),
            "patient_ids": patient_ids,
            "sample_ids": sample_ids,
            "patient_sample_map": [
                {"patient_id": patient_id, "sample_id": sample_id}
                for patient_id, sample_id in zip(patient_ids, sample_ids, strict=True)
            ],
            "pam50_counts": {label: counts.get(label, 0) for label in PAM50_CLASS_ORDER},
        }
    return partitions


def _validation_results(
    records: Sequence[Mapping[str, str]],
    labels: Sequence[str],
    partitions: Mapping[str, Mapping[str, Any]],
    target_sizes: Mapping[str, int],
    matrix_header_validation: Mapping[str, Any],
) -> dict[str, Any]:
    cohort_patients = [row["patient_id"] for row in records]
    cohort_samples = [row["sample_id"] for row in records]
    cohort_label_counts = Counter(labels)
    patient_occurrences = Counter(
        patient_id
        for partition in PARTITION_ORDER
        for patient_id in partitions[partition]["patient_ids"]
    )
    sample_occurrences = Counter(
        sample_id
        for partition in PARTITION_ORDER
        for sample_id in partitions[partition]["sample_ids"]
    )
    disjoint_patients = all(
        set(partitions[left]["patient_ids"]).isdisjoint(partitions[right]["patient_ids"])
        for index, left in enumerate(PARTITION_ORDER)
        for right in PARTITION_ORDER[index + 1 :]
    )
    disjoint_samples = all(
        set(partitions[left]["sample_ids"]).isdisjoint(partitions[right]["sample_ids"])
        for index, left in enumerate(PARTITION_ORDER)
        for right in PARTITION_ORDER[index + 1 :]
    )
    patient_partition_exact = set(patient_occurrences) == set(cohort_patients) and all(
        count == 1 for count in patient_occurrences.values()
    )
    sample_partition_exact = set(sample_occurrences) == set(cohort_samples) and all(
        count == 1 for count in sample_occurrences.values()
    )
    class_counts_reconcile = all(
        sum(partitions[partition]["pam50_counts"][label] for partition in PARTITION_ORDER)
        == cohort_label_counts[label]
        for label in PAM50_CLASS_ORDER
    )
    target_sizes_match = all(
        partitions[partition]["n_patients"] == target_sizes[partition]
        for partition in PARTITION_ORDER
    )
    stratified_class_coverage = all(
        partitions[partition]["pam50_counts"][label] >= 1
        for partition in PARTITION_ORDER
        for label in PAM50_CLASS_ORDER
    )

    checks = [
        {
            "name": "unique_patient_ids",
            "passed": len(cohort_patients) == len(set(cohort_patients)),
            "observed": len(set(cohort_patients)),
        },
        {
            "name": "unique_sample_ids",
            "passed": len(cohort_samples) == len(set(cohort_samples)),
            "observed": len(set(cohort_samples)),
        },
        {
            "name": "all_cohort_samples_match_expression_matrix_header",
            "passed": bool(matrix_header_validation["sample_order_matches_cohort"]),
            "observed": matrix_header_validation["sample_column_count"],
        },
        {
            "name": "patient_sets_pairwise_disjoint",
            "passed": disjoint_patients,
        },
        {
            "name": "sample_sets_pairwise_disjoint",
            "passed": disjoint_samples,
        },
        {
            "name": "every_cohort_patient_appears_exactly_once",
            "passed": patient_partition_exact,
            "observed": sum(patient_occurrences.values()),
        },
        {
            "name": "every_cohort_sample_appears_exactly_once",
            "passed": sample_partition_exact,
            "observed": sum(sample_occurrences.values()),
        },
        {
            "name": "pam50_counts_reconcile_with_cohort",
            "passed": class_counts_reconcile,
        },
        {
            "name": "largest_remainder_partition_sizes_match",
            "passed": target_sizes_match,
            "target_sizes": dict(target_sizes),
        },
        {
            "name": "every_pam50_class_represented_in_each_partition",
            "passed": stratified_class_coverage,
        },
    ]
    return {
        "passed": all(check["passed"] for check in checks),
        "checks": checks,
        "cohort_patient_count": len(cohort_patients),
        "cohort_sample_count": len(cohort_samples),
        "cohort_pam50_counts": {label: cohort_label_counts[label] for label in PAM50_CLASS_ORDER},
    }


def generate_split_manifest(
    *,
    cohort_path: Path,
    expression_matrix_path: Path,
    matrix_manifest_path: Path,
    config_path: Path,
    project_root: Path,
    expected_patient_count: int,
    expected_gene_count: int,
    expression_scale: str,
    random_seed: int,
    proportions: Mapping[str, float],
    rounding_tie_order: Sequence[str],
    test_random_state_offset: int = 0,
    validation_random_state_offset: int = 1,
) -> dict[str, Any]:
    """Create a deterministic split manifest without consulting expression values."""
    project_root = Path(project_root).resolve()
    cohort_path = Path(cohort_path).resolve()
    expression_matrix_path = Path(expression_matrix_path).resolve()
    matrix_manifest_path = Path(matrix_manifest_path).resolve()
    config_path = Path(config_path).resolve()
    if expected_patient_count < 1 or expected_gene_count < 1:
        raise ValueError("Expected patient and gene counts must be positive")
    if not 0 <= random_seed < 2**32:
        raise ValueError("Random seed must be an integer in [0, 2**32)")
    if expression_scale.strip() != "log2(normalized_count + 1)":
        raise ValueError("Phase 1F requires expression units log2(normalized_count + 1)")

    records, labels = _read_validated_cohort(
        cohort_path,
        expected_patient_count=expected_patient_count,
    )
    if len(records) != expected_patient_count:
        raise ValueError(
            f"Expected {expected_patient_count} approved patients, found {len(records)}"
        )
    sample_ids = [row["sample_id"] for row in records]
    matrix_header_validation = _read_expression_matrix_header(
        expression_matrix_path,
        matrix_manifest_path,
        project_root=project_root,
        cohort_sample_ids=sample_ids,
        expected_gene_count=expected_gene_count,
        expression_scale=expression_scale,
    )

    target_sizes = largest_remainder_sizes(
        len(records),
        proportions,
        rounding_tie_order,
    )
    indices_by_partition = _stratified_partition_indices(
        labels,
        target_sizes=target_sizes,
        random_seed=random_seed,
        test_seed_offset=test_random_state_offset,
        validation_seed_offset=validation_random_state_offset,
    )
    partitions = _partition_records(records, labels, indices_by_partition)
    validation = _validation_results(
        records,
        labels,
        partitions,
        target_sizes,
        matrix_header_validation,
    )
    if not validation["passed"]:
        failed = [check["name"] for check in validation["checks"] if not check["passed"]]
        raise ValueError(f"Generated split failed validation checks: {failed}")

    cohort_counts = Counter(labels)
    matrix_manifest = json.loads(matrix_manifest_path.read_text(encoding="utf-8"))
    return {
        "schema_version": 1,
        "split_id": "split_v1",
        "created_at": datetime.now(ZoneInfo("Asia/Singapore")).isoformat(timespec="seconds"),
        "random_seed": random_seed,
        "split_proportions": {key: float(proportions[key]) for key in PARTITION_ORDER},
        "split_unit": "patient_id",
        "stratification_label": "pam50_normalized_label",
        "cohort_pam50_counts": {label: cohort_counts[label] for label in PAM50_CLASS_ORDER},
        "target_patient_counts": target_sizes,
        "rounding": {
            "method": "largest_remainder",
            "description": (
                "Floor n*p for each partition, then assign remaining patients by "
                "descending fractional remainder; ties use rounding_tie_order."
            ),
            "tie_break_order": list(rounding_tie_order),
        },
        "generation_method": {
            "algorithm": "two-stage sklearn.model_selection.train_test_split",
            "stratified": True,
            "test_random_state": random_seed + test_random_state_offset,
            "validation_random_state": random_seed + validation_random_state_offset,
            "test_random_state_offset": test_random_state_offset,
            "validation_random_state_offset": validation_random_state_offset,
            "library": "scikit-learn",
            "library_version": sklearn.__version__,
            "partition_sequence": [
                "stratified test set from all approved patients",
                "stratified validation set from the non-test patients",
                "training set is the remaining patients",
            ],
            "serialization_order": (
                "members are serialized in original cohort row order after assignment; "
                "this does not change partition membership"
            ),
            "expression_values_used_for_assignment": False,
        },
        "source_cohort": {
            "path": _manifest_path(cohort_path, project_root),
            "sha256": sha256_file(cohort_path),
            "patient_count": len(records),
            "sample_count": len(sample_ids),
        },
        "expression_matrix_header_validation": {
            "path": _manifest_path(expression_matrix_path, project_root),
            "matrix_manifest_path": _manifest_path(matrix_manifest_path, project_root),
            "matrix_manifest_sha256": sha256_file(matrix_manifest_path),
            "matrix_sha256_from_manifest": matrix_manifest["matrix"].get("sha256"),
            **matrix_header_validation,
        },
        "configuration": {
            "path": _manifest_path(config_path, project_root),
            "sha256": sha256_file(config_path),
            "random_seed": random_seed,
            "split_proportions": {key: float(proportions[key]) for key in PARTITION_ORDER},
            "expected_patient_count": expected_patient_count,
            "expected_gene_count": expected_gene_count,
        },
        "partitions": partitions,
        "validation_results": validation,
        "test_set_policy": (
            "locked during model development; do not inspect test-set performance or use it "
            "for model, feature, or hyperparameter selection"
        ),
    }


def write_frozen_split_manifest(manifest: Mapping[str, Any], output_path: Path) -> bool:
    """Write split_v1 once; identical reruns validate but never rewrite the frozen file.

    Returns True if newly created, False if an identical frozen manifest already exists.
    A changed cohort, seed, configuration, library version, or partition assignment requires
    a new split version instead of overwriting split_v1.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if output_path.exists():
        try:
            existing = json.loads(output_path.read_text(encoding="utf-8"))
        except json.JSONDecodeError as exc:
            raise ValueError(
                f"Existing frozen split manifest is invalid JSON: {output_path}"
            ) from exc
        existing_without_time = dict(existing)
        generated_without_time = dict(manifest)
        existing_without_time.pop("created_at", None)
        generated_without_time.pop("created_at", None)
        if existing_without_time == generated_without_time:
            return False
        raise FileExistsError(
            f"Refusing to overwrite frozen split manifest {output_path}; "
            "create a new split version for changed inputs or assignments"
        )

    temp_path = output_path.with_name(output_path.name + ".tmp")
    try:
        with temp_path.open("w", encoding="utf-8", newline="\n") as handle:
            json.dump(manifest, handle, indent=2, sort_keys=True)
            handle.write("\n")
        temp_path.replace(output_path)
    finally:
        if temp_path.exists():
            temp_path.unlink()
    return True
