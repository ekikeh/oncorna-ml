"""Read-only split checks and a fold-fitted expression-prevalence filter.

Input expression is already log2(normalized_count + 1). This module never
normalizes it or uses labels to decide which genes to retain.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np
import pandas as pd
import yaml
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.utils.validation import check_is_fitted

from oncorna.qc import MISSING_TOKENS


@dataclass(frozen=True)
class AllGeneView:
    """Unfiltered cohort values in samples-by-source-genes order."""

    gene_ids: tuple[str, ...]
    sample_ids: tuple[str, ...]
    values: np.ndarray
    source_sample_count: int


def read_all_gene_view(
    path: Path,
    cohort_sample_ids: Sequence[str],
    *,
    expected_gene_count: int,
) -> AllGeneView:
    """Copy numeric source values for approved samples without gene filtering.

    Float64 preserves the parsed source precision around the strict >1 boundary.
    The source's existing log2-normalized values receive no new transformation.
    """
    sample_ids = tuple(cohort_sample_ids)
    if not sample_ids or len(set(sample_ids)) != len(sample_ids):
        raise ValueError("Approved sample IDs must be nonempty and unique")
    if expected_gene_count < 1:
        raise ValueError("Expected source gene count must be positive")
    values = np.empty((len(sample_ids), expected_gene_count), dtype=np.float64)
    gene_ids: list[str] = []
    seen_genes: set[str] = set()
    with gzip.open(path, "rt", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, None)
        if not header or header[0] != "sample" or len(header) != len(set(header)):
            raise ValueError("Source header must start with sample and have unique columns")
        source_samples = header[1:]
        index_by_sample = {sample: index + 1 for index, sample in enumerate(source_samples)}
        missing = set(sample_ids) - set(index_by_sample)
        if missing:
            raise ValueError(
                f"Approved samples missing from expression source: {sorted(missing)[:5]}"
            )
        selected_indices = [index_by_sample[sample] for sample in sample_ids]
        for gene_index, row in enumerate(reader):
            if gene_index >= expected_gene_count or len(row) != len(header):
                raise ValueError(
                    "Expression source gene count or row width differs from expectation"
                )
            gene_id = row[0]
            if not gene_id or gene_id != gene_id.strip() or gene_id in seen_genes:
                raise ValueError(f"Invalid or duplicate source gene ID: {gene_id!r}")
            seen_genes.add(gene_id)
            gene_ids.append(gene_id)
            for sample_index, source_index in enumerate(selected_indices):
                token = row[source_index].strip()
                if token.casefold() in MISSING_TOKENS:
                    value = math.nan
                else:
                    try:
                        value = float(token)
                    except ValueError as exc:
                        raise ValueError(
                            f"Nonnumeric expression for gene {gene_id}, "
                            f"sample {sample_ids[sample_index]}"
                        ) from exc
                if math.isinf(value):
                    raise ValueError(
                        f"Infinite expression for gene {gene_id}, sample {sample_ids[sample_index]}"
                    )
                values[sample_index, gene_index] = value
    if len(gene_ids) != expected_gene_count:
        raise ValueError(f"Expected {expected_gene_count} source genes, found {len(gene_ids)}")
    return AllGeneView(tuple(gene_ids), sample_ids, values, len(source_samples))


def partition_frames(
    view: AllGeneView,
    partitions: Mapping[str, Sequence[tuple[str, str]]],
) -> dict[str, pd.DataFrame]:
    """Select frozen patient/sample rows without computing expression statistics."""
    if view.values.shape != (len(view.sample_ids), len(view.gene_ids)):
        raise ValueError("All-gene view dimensions disagree with its identifiers")
    source_positions = {sample: index for index, sample in enumerate(view.sample_ids)}
    if len(source_positions) != len(view.sample_ids):
        raise ValueError("All-gene view has duplicate sample IDs")
    result: dict[str, pd.DataFrame] = {}
    for name in ("train", "validation", "test"):
        sample_ids = [sample for _, sample in partitions[name]]
        if any(sample not in source_positions for sample in sample_ids):
            raise ValueError(f"{name} sample missing from the all-gene view")
        row_indices = [source_positions[sample] for sample in sample_ids]
        result[name] = pd.DataFrame(
            view.values[row_indices, :], index=sample_ids, columns=view.gene_ids
        )
    return result


def sha256_file(path: Path) -> str:
    """Hash a file's actual bytes, including its line endings."""
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def identity_digest(pairs: Sequence[tuple[str, str]]) -> str:
    """Hash ordered patient/sample pairs with a documented UTF-8 encoding."""
    payload = json.dumps(list(pairs), ensure_ascii=False, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def verify_historical_config(path: Path, recorded_sha256: str) -> dict[str, str]:
    """Accept only exact bytes or a provable LF/CRLF-only variant.

    The recorded raw digest is never replaced. The returned mode makes any
    Windows checkout difference visible to the caller and its manifest.
    """
    raw = Path(path).read_bytes()
    actual = hashlib.sha256(raw).hexdigest()
    if actual == recorded_sha256:
        mode = "exact_raw_bytes"
    else:
        lf = raw.replace(b"\r\n", b"\n")
        crlf = lf.replace(b"\n", b"\r\n")
        variants = {
            "recorded_lf_current_crlf": hashlib.sha256(lf).hexdigest(),
            "recorded_crlf_current_lf": hashlib.sha256(crlf).hexdigest(),
        }
        matching = [name for name, digest in variants.items() if digest == recorded_sha256]
        if len(matching) != 1 or b"\r" in lf:
            raise ValueError("Frozen split configuration differs beyond LF/CRLF line endings")
        mode = matching[0]
    return {
        "path": str(path),
        "recorded_raw_sha256": recorded_sha256,
        "current_raw_sha256": actual,
        "current_lf_sha256": hashlib.sha256(raw.replace(b"\r\n", b"\n")).hexdigest(),
        "compatibility_mode": mode,
    }


def validate_frozen_partitions(
    cohort_records: Sequence[Mapping[str, str]],
    split: Mapping[str, Any],
) -> dict[str, tuple[tuple[str, str], ...]]:
    """Validate and return the frozen partitions in recorded cohort row order."""
    if split.get("schema_version") != 1 or split.get("split_id") != "split_v1":
        raise ValueError("Expected the existing split_v1 schema")
    if split.get("split_unit") != "patient_id":
        raise ValueError("Frozen split must assign whole patients")
    if split.get("validation_results", {}).get("passed") is not True:
        raise ValueError("Frozen split does not record passing validation")

    cohort_pairs = tuple((row["patient_id"], row["sample_id"]) for row in cohort_records)
    if not cohort_pairs or len(set(cohort_pairs)) != len(cohort_pairs):
        raise ValueError("Approved cohort has missing or duplicate patient/sample pairs")
    if len({patient for patient, _ in cohort_pairs}) != len(cohort_pairs):
        raise ValueError("Approved cohort must have one sample per patient")
    if len({sample for _, sample in cohort_pairs}) != len(cohort_pairs):
        raise ValueError("Approved cohort has duplicate sample IDs")
    position = {pair: index for index, pair in enumerate(cohort_pairs)}
    labels_by_pair = {
        (row["patient_id"], row["sample_id"]): row["pam50_normalized_label"]
        for row in cohort_records
    }
    source = split.get("source_cohort", {})
    if source.get("patient_count") != len(cohort_pairs) or source.get("sample_count") != len(
        cohort_pairs
    ):
        raise ValueError("Frozen split cohort counts differ from the approved cohort")

    result: dict[str, tuple[tuple[str, str], ...]] = {}
    seen: set[tuple[str, str]] = set()
    for name in ("train", "validation", "test"):
        part = split.get("partitions", {}).get(name)
        if not isinstance(part, dict):
            raise ValueError(f"Frozen split is missing {name}")
        pairs = tuple(
            (item["patient_id"], item["sample_id"]) for item in part["patient_sample_map"]
        )
        if tuple(part["patient_ids"]) != tuple(patient for patient, _ in pairs):
            raise ValueError(f"{name} patient IDs disagree with its patient/sample map")
        if tuple(part["sample_ids"]) != tuple(sample for _, sample in pairs):
            raise ValueError(f"{name} sample IDs disagree with its patient/sample map")
        if len(pairs) != part["n_patients"] or len(pairs) != part["n_samples"]:
            raise ValueError(f"{name} partition count is inconsistent")
        if any(pair not in position for pair in pairs):
            raise ValueError(f"{name} contains a patient/sample pair outside the cohort")
        positions = [position[pair] for pair in pairs]
        if positions != sorted(positions):
            raise ValueError(f"{name} order differs from the approved cohort order")
        if seen.intersection(pairs):
            raise ValueError("Frozen patient/sample partitions overlap")
        if dict(Counter(labels_by_pair[pair] for pair in pairs)) != part["pam50_counts"]:
            raise ValueError(f"{name} subtype counts differ from the approved cohort")
        seen.update(pairs)
        result[name] = pairs
    if seen != set(cohort_pairs):
        raise ValueError("Frozen partitions do not cover the approved cohort exactly")
    if dict(Counter(labels_by_pair.values())) != split.get("cohort_pam50_counts"):
        raise ValueError("Frozen cohort subtype counts differ from the approved cohort")
    return result


def verify_split_provenance(
    split: Mapping[str, Any],
    *,
    cohort_path: Path,
    historical_config_path: Path,
    historical_matrix_path: Path,
    historical_matrix_manifest_path: Path,
) -> dict[str, Any]:
    """Check original input bytes and explain a config line-ending difference."""
    cohort_sha256 = sha256_file(cohort_path)
    if cohort_sha256 != split["source_cohort"]["sha256"]:
        raise ValueError("Approved cohort checksum differs from the frozen split")
    config = verify_historical_config(historical_config_path, split["configuration"]["sha256"])
    with Path(historical_config_path).open("r", encoding="utf-8") as handle:
        parsed = yaml.safe_load(handle)
    if (
        parsed["random_seed"] != split["random_seed"]
        or parsed["splitting"]["split_unit"] != "patient_id"
        or parsed["splitting"]["stratification_label"] != "pam50_normalized_label"
        or parsed["splitting"]["proportions"] != split["split_proportions"]
        or parsed["qc"]["expected_cohort_samples"]
        != split["configuration"]["expected_patient_count"]
        or parsed["ml_matrix"]["expected_gene_count"]
        != split["configuration"]["expected_gene_count"]
    ):
        raise ValueError("Parsed project configuration disagrees with frozen split settings")

    matrix_record = split["expression_matrix_header_validation"]
    matrix_manifest_sha256 = sha256_file(historical_matrix_manifest_path)
    if matrix_manifest_sha256 != matrix_record["matrix_manifest_sha256"]:
        raise ValueError("Historical matrix manifest checksum differs from the frozen split")
    with Path(historical_matrix_manifest_path).open("r", encoding="utf-8") as handle:
        matrix_manifest = json.load(handle)
    matrix_sha256 = sha256_file(historical_matrix_path)
    if (
        matrix_sha256 != matrix_record["matrix_sha256_from_manifest"]
        or matrix_sha256 != matrix_manifest["matrix"]["sha256"]
    ):
        raise ValueError("Historical matrix checksum differs from its recorded provenance")
    return {
        "cohort_sha256": cohort_sha256,
        "historical_config": config,
        "historical_matrix_sha256": matrix_sha256,
        "historical_matrix_manifest_sha256": matrix_manifest_sha256,
    }


class GenePrevalenceFilter(TransformerMixin, BaseEstimator):
    """Keep genes above a fixed expression threshold in enough fitting samples.

    X must be a numeric pandas DataFrame with samples as rows and unique gene
    identifiers as columns. ``fit`` ignores y; ``transform`` reuses its saved
    mask and requires the full input gene schema in the same order.
    """

    def __init__(self, expression_gt: float = 1.0, minimum_sample_fraction: float = 0.20):
        self.expression_gt = expression_gt
        self.minimum_sample_fraction = minimum_sample_fraction

    @staticmethod
    def _validated_values(X: pd.DataFrame) -> np.ndarray:
        if not isinstance(X, pd.DataFrame):
            raise TypeError("GenePrevalenceFilter requires a pandas DataFrame with gene IDs")
        if X.empty or X.columns.has_duplicates or X.index.has_duplicates:
            raise ValueError("Expression input must be nonempty with unique gene and sample IDs")
        if any(not isinstance(gene, str) or not gene or gene != gene.strip() for gene in X.columns):
            raise ValueError("Every input gene ID must be a nonblank, unpadded string")
        if not all(pd.api.types.is_numeric_dtype(dtype) for dtype in X.dtypes):
            raise ValueError("Expression input must contain only numeric gene columns")
        values = X.to_numpy(dtype=np.float64, na_value=np.nan, copy=False)
        if np.isinf(values).any():
            raise ValueError("Infinite expression values are not permitted")
        return values

    def fit(self, X: pd.DataFrame, y: Any = None) -> GenePrevalenceFilter:
        """Learn a label-independent mask from this call's samples only."""
        if not math.isfinite(self.expression_gt):
            raise ValueError("expression_gt must be finite")
        if not math.isfinite(self.minimum_sample_fraction) or not (
            0 < self.minimum_sample_fraction <= 1
        ):
            raise ValueError("minimum_sample_fraction must be in (0, 1]")
        values = self._validated_values(X)
        required = max(1, math.ceil(self.minimum_sample_fraction * len(X)))
        counts = np.count_nonzero(np.isfinite(values) & (values > self.expression_gt), axis=0)
        self.feature_names_in_ = np.asarray(X.columns, dtype=object)
        self.n_features_in_ = X.shape[1]
        self.n_samples_fit_ = len(X)
        self.required_count_ = required
        self.prevalence_counts_ = counts
        self.support_mask_ = counts >= required
        self.selected_gene_ids_ = tuple(self.feature_names_in_[self.support_mask_])
        self.fitting_sample_ids_ = tuple(str(sample_id) for sample_id in X.index)
        self.missing_value_count_ = int(np.isnan(values).sum())
        self.infinite_value_count_ = 0
        return self

    def transform(self, X: pd.DataFrame) -> pd.DataFrame:
        """Select the fitted genes without learning from these samples."""
        check_is_fitted(self, "support_mask_")
        self._validated_values(X)
        if tuple(X.columns) != tuple(self.feature_names_in_):
            raise ValueError("Incoming gene identifiers or order differ from the fitted schema")
        return X.loc[:, list(self.selected_gene_ids_)].copy()

    def get_support(self, indices: bool = False) -> np.ndarray:
        check_is_fitted(self, "support_mask_")
        return np.flatnonzero(self.support_mask_) if indices else self.support_mask_.copy()

    def get_feature_names_out(self, input_features: Sequence[str] | None = None) -> np.ndarray:
        check_is_fitted(self, "support_mask_")
        if input_features is not None and tuple(input_features) != tuple(self.feature_names_in_):
            raise ValueError("Input feature names differ from the fitted schema")
        return np.asarray(self.selected_gene_ids_, dtype=object)
