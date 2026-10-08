"""Build a reproducible, metadata-free ML expression matrix from Phase 1C/1D outputs.

Source values are copied verbatim from the HiSeqV2 file. They remain
log2(normalized_count + 1); this module does not log-transform, normalize,
scale, impute, or use PAM50 labels to select genes.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import io
import json
import math
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any, Mapping, Sequence
from zoneinfo import ZoneInfo

from oncorna.qc import read_classification_cohort

GENE_LIST_FIELDS = (
    "gene_id",
    "selected_filter_expression_gt",
    "n_samples_above_selected_threshold",
    "fraction_samples_above_selected_threshold",
    "required_samples_for_selected_filter",
)

SAMPLE_METRIC_REQUIRED_FIELDS = (
    "sample_id",
    "sample_included_in_qc_cohort",
    "sample_exclusion_reason",
    "tukey_outlier_aggregate_signal",
    "n_missing_values",
    "n_infinite_values",
)

GENE_SUMMARY_REQUIRED_FIELDS = (
    "gene_id",
    "n_samples_in_qc_cohort",
    "n_samples_above_selected_threshold",
    "fraction_samples_above_selected_threshold",
    "selected_filter_expression_gt",
    "required_samples_for_selected_filter",
    "retained_by_prevalence_filter",
)

MISSING_TOKENS = {"", "na", "nan", "null", "none", "n/a", "--"}


@dataclass(frozen=True)
class MLMatrixBuildResult:
    """Counts and output paths from a Phase 1E matrix build."""

    gene_ids: tuple[str, ...]
    sample_ids: tuple[str, ...]
    source_gene_count: int
    source_sample_count: int
    noncohort_source_sample_count: int
    required_sample_count: int
    flagged_outlier_sample_ids: tuple[str, ...]
    matrix_path: Path
    gene_list_path: Path
    manifest_path: Path
    matrix_sha256: str


def sha256_file(path: Path, *, chunk_size: int = 1024 * 1024) -> str:
    """Return a streaming SHA-256 digest for a file."""
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _read_tsv(path: Path, required_fields: Sequence[str]) -> list[dict[str, str]]:
    path = Path(path)
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        actual_fields = set(reader.fieldnames or ())
        missing = sorted(set(required_fields) - actual_fields)
        if missing:
            raise ValueError(f"{path} is missing required columns: {missing}")
        return [
            {key: (value or "").strip() for key, value in row.items() if key is not None}
            for row in reader
        ]


def _as_bool(value: str, *, field: str, row_id: str) -> bool:
    normalized = value.strip().casefold()
    if normalized in {"true", "1"}:
        return True
    if normalized in {"false", "0"}:
        return False
    raise ValueError(f"Invalid boolean {value!r} in {field} for {row_id}")


def _as_int(value: str, *, field: str, row_id: str) -> int:
    try:
        return int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Invalid integer {value!r} in {field} for {row_id}") from exc


def _as_float(value: str, *, field: str, row_id: str) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Invalid number {value!r} in {field} for {row_id}") from exc
    if not math.isfinite(result):
        raise ValueError(f"Non-finite value {value!r} in {field} for {row_id}")
    return result


def _load_sample_qc(path: Path, sample_ids: Sequence[str]) -> list[str]:
    """Validate the Phase 1D sample dispositions without using PAM50 labels."""
    rows = _read_tsv(path, SAMPLE_METRIC_REQUIRED_FIELDS)
    by_id: dict[str, dict[str, str]] = {}
    for row in rows:
        sample_id = row["sample_id"]
        if not sample_id:
            raise ValueError("Phase 1D sample metrics contain a blank sample_id")
        if sample_id in by_id:
            raise ValueError(f"Duplicate sample_id in Phase 1D sample metrics: {sample_id}")
        by_id[sample_id] = row
    if set(by_id) != set(sample_ids):
        missing = sorted(set(sample_ids) - set(by_id))
        extra = sorted(set(by_id) - set(sample_ids))
        raise ValueError(
            "Phase 1D sample metrics do not match the classification cohort: "
            f"missing={missing[:5]}, extra={extra[:5]}"
        )

    flagged_outliers: list[str] = []
    for sample_id in sample_ids:
        row = by_id[sample_id]
        included = _as_bool(
            row["sample_included_in_qc_cohort"],
            field="sample_included_in_qc_cohort",
            row_id=sample_id,
        )
        if not included or row["sample_exclusion_reason"] != "retained_for_qc":
            raise ValueError(f"Phase 1C cohort sample did not pass Phase 1D sample QC: {sample_id}")
        missing_count = _as_int(row["n_missing_values"], field="n_missing_values", row_id=sample_id)
        infinite_count = _as_int(
            row["n_infinite_values"], field="n_infinite_values", row_id=sample_id
        )
        if missing_count or infinite_count:
            raise ValueError(
                f"Phase 1D reports non-finite expression values for retained sample {sample_id}: "
                f"missing={missing_count}, infinite={infinite_count}"
            )
        if _as_bool(
            row["tukey_outlier_aggregate_signal"],
            field="tukey_outlier_aggregate_signal",
            row_id=sample_id,
        ):
            flagged_outliers.append(sample_id)
    return flagged_outliers


def _load_gene_filter_results(
    path: Path,
    *,
    expected_sample_count: int,
    expected_source_gene_count: int | None,
    expression_threshold: float,
    minimum_sample_fraction: float,
    expected_gene_count: int | None,
) -> tuple[list[dict[str, str]], list[dict[str, Any]], int]:
    """Validate the Phase 1D prevalence flags against the configured rule."""
    rows = _read_tsv(path, GENE_SUMMARY_REQUIRED_FIELDS)
    if not rows:
        raise ValueError("Phase 1D gene summary contains no genes")
    if expected_source_gene_count is not None and len(rows) != expected_source_gene_count:
        raise ValueError(
            "Phase 1D gene-summary row count does not match configuration: "
            f"expected {expected_source_gene_count}, found {len(rows)}"
        )

    gene_ids: set[str] = set()
    required_count = max(1, math.ceil(minimum_sample_fraction * expected_sample_count))
    selected_rows: list[dict[str, Any]] = []
    for row in rows:
        gene_id = row["gene_id"]
        if not gene_id:
            raise ValueError("Phase 1D gene summary contains a blank gene_id")
        if gene_id in gene_ids:
            raise ValueError(f"Duplicate gene_id in Phase 1D gene summary: {gene_id}")
        gene_ids.add(gene_id)
        n_samples = _as_int(
            row["n_samples_in_qc_cohort"], field="n_samples_in_qc_cohort", row_id=gene_id
        )
        n_above = _as_int(
            row["n_samples_above_selected_threshold"],
            field="n_samples_above_selected_threshold",
            row_id=gene_id,
        )
        stored_required = _as_int(
            row["required_samples_for_selected_filter"],
            field="required_samples_for_selected_filter",
            row_id=gene_id,
        )
        stored_threshold = _as_float(
            row["selected_filter_expression_gt"],
            field="selected_filter_expression_gt",
            row_id=gene_id,
        )
        fraction = _as_float(
            row["fraction_samples_above_selected_threshold"],
            field="fraction_samples_above_selected_threshold",
            row_id=gene_id,
        )
        stored_keep = _as_bool(
            row["retained_by_prevalence_filter"],
            field="retained_by_prevalence_filter",
            row_id=gene_id,
        )
        if n_samples != expected_sample_count:
            raise ValueError(
                f"Gene {gene_id} was summarized over {n_samples} samples; "
                f"expected {expected_sample_count}"
            )
        if stored_required != required_count:
            raise ValueError(
                f"Gene {gene_id} records {stored_required} required samples; "
                f"configured rule requires {required_count}"
            )
        if not math.isclose(stored_threshold, expression_threshold, rel_tol=0, abs_tol=1e-10):
            raise ValueError(
                f"Gene {gene_id} filter threshold {stored_threshold} does not match "
                f"configured threshold {expression_threshold}"
            )
        if not 0 <= n_above <= expected_sample_count:
            raise ValueError(f"Invalid prevalence count for gene {gene_id}: {n_above}")
        if not math.isclose(
            fraction,
            n_above / expected_sample_count,
            rel_tol=0,
            abs_tol=1e-8,
        ):
            raise ValueError(f"Prevalence fraction disagrees with count for gene {gene_id}")
        calculated_keep = n_above >= required_count
        if stored_keep != calculated_keep:
            raise ValueError(
                f"Prevalence flag for {gene_id} disagrees with the configured threshold"
            )
        if calculated_keep:
            selected_rows.append(
                {
                    "gene_id": gene_id,
                    "selected_filter_expression_gt": stored_threshold,
                    "n_samples_above_selected_threshold": n_above,
                    "fraction_samples_above_selected_threshold": fraction,
                    "required_samples_for_selected_filter": required_count,
                }
            )

    if expected_gene_count is not None and len(selected_rows) != expected_gene_count:
        raise ValueError(
            "Retained-gene count does not match the Phase 1E configuration: "
            f"expected {expected_gene_count}, found {len(selected_rows)}"
        )
    return rows, selected_rows, required_count


def _write_tsv(path: Path, rows: Sequence[Mapping[str, Any]], fields: Sequence[str]) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp_path = path.with_name(path.name + ".tmp")
    try:
        with temp_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=list(fields),
                delimiter="\t",
                extrasaction="ignore",
                lineterminator="\n",
            )
            writer.writeheader()
            writer.writerows(rows)
        temp_path.replace(path)
    finally:
        if temp_path.exists():
            temp_path.unlink()


def _verify_output_dimensions(
    path: Path,
    *,
    sample_ids: Sequence[str],
    gene_ids: Sequence[str],
) -> None:
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, None)
        expected_header = ["gene_id", *sample_ids]
        if header != expected_header:
            raise ValueError("Generated ML matrix header does not match the manifest sample order")
        seen: set[str] = set()
        observed: list[str] = []
        for row_number, row in enumerate(reader, start=2):
            if len(row) != len(expected_header):
                raise ValueError(
                    f"Generated ML matrix row {row_number} has {len(row)} fields; "
                    f"expected {len(expected_header)}"
                )
            gene_id = row[0]
            if not gene_id or gene_id in seen:
                raise ValueError(
                    f"Generated ML matrix has a blank or duplicate gene ID: {gene_id!r}"
                )
            seen.add(gene_id)
            observed.append(gene_id)
        if observed != list(gene_ids):
            raise ValueError("Generated ML matrix gene order does not match selected gene list")


def _manifest_path(path: Path, project_root: Path) -> str:
    try:
        return path.relative_to(project_root).as_posix()
    except ValueError:
        return path.as_posix()


def _input_record(path: Path, project_root: Path) -> dict[str, str]:
    path = Path(path).resolve()
    return {
        "path": _manifest_path(path, project_root),
        "sha256": sha256_file(path),
    }


def build_ml_matrix(
    *,
    expression_path: Path,
    cohort_path: Path,
    gene_summary_path: Path,
    sample_metrics_path: Path,
    config_path: Path,
    output_path: Path,
    gene_list_path: Path,
    manifest_path: Path,
    project_root: Path,
    expression_scale: str,
    expression_threshold: float,
    minimum_sample_fraction: float,
    expected_sample_count: int,
    expected_gene_count: int,
    expected_source_gene_count: int | None,
    expression_transform: str,
    config_snapshot: Mapping[str, Any],
) -> MLMatrixBuildResult:
    """Stream the approved sample/gene subset into a verbatim-value matrix."""
    project_root = Path(project_root).resolve()
    expression_path = Path(expression_path).resolve()
    cohort_path = Path(cohort_path).resolve()
    gene_summary_path = Path(gene_summary_path).resolve()
    sample_metrics_path = Path(sample_metrics_path).resolve()
    config_path = Path(config_path).resolve()
    output_path = Path(output_path).resolve()
    gene_list_path = Path(gene_list_path).resolve()
    manifest_path = Path(manifest_path).resolve()

    if expected_sample_count < 1 or expected_gene_count < 1:
        raise ValueError("Expected sample and retained-gene counts must be positive")
    if expected_source_gene_count is not None and expected_source_gene_count < 1:
        raise ValueError("Expected source-gene count must be positive")
    if expected_source_gene_count is not None and expected_gene_count > expected_source_gene_count:
        raise ValueError("Retained-gene count cannot exceed source-gene count")
    if not math.isfinite(expression_threshold):
        raise ValueError("Expression prevalence threshold must be finite")
    if not math.isfinite(minimum_sample_fraction) or not 0 < minimum_sample_fraction <= 1:
        raise ValueError("Minimum sample fraction must be greater than 0 and at most 1")
    if expression_scale.strip() != "log2(normalized_count + 1)":
        raise ValueError("Phase 1E requires expression units log2(normalized_count + 1)")

    processed_dir = (project_root / "data" / "processed").resolve()
    output_paths = (output_path, gene_list_path, manifest_path)
    if any(not path.is_relative_to(processed_dir) for path in output_paths):
        raise ValueError("Phase 1E outputs must stay under Git-ignored data/processed/")

    if expression_transform.strip().casefold() not in {
        "none",
        "none; copy selected source values verbatim",
    }:
        raise ValueError(
            "Phase 1E does not apply a second expression transformation; "
            f"configured expression_transform was {expression_transform!r}"
        )
    cohort_records = read_classification_cohort(
        cohort_path,
        expected_sample_count=expected_sample_count,
    )
    sample_ids = [row["sample_id"] for row in cohort_records]
    if len(sample_ids) != len(set(sample_ids)):
        raise ValueError("Classification cohort has duplicate sample IDs")
    patient_ids = [row["patient_id"] for row in cohort_records]
    if len(patient_ids) != len(set(patient_ids)):
        raise ValueError("Classification cohort has duplicate patient IDs")

    flagged_outliers = _load_sample_qc(sample_metrics_path, sample_ids)
    all_gene_rows, selected_gene_rows, required_samples = _load_gene_filter_results(
        gene_summary_path,
        expected_sample_count=expected_sample_count,
        expected_source_gene_count=expected_source_gene_count,
        expression_threshold=expression_threshold,
        minimum_sample_fraction=minimum_sample_fraction,
        expected_gene_count=expected_gene_count,
    )
    selected_gene_ids = [row["gene_id"] for row in selected_gene_rows]

    output_path.parent.mkdir(parents=True, exist_ok=True)
    temp_matrix_path = output_path.with_name(output_path.name + ".tmp")
    source_gene_count = 0
    source_sample_count = 0
    noncohort_source_sample_count = 0
    seen_source_genes: set[str] = set()
    selected_index = 0
    sample_indices: list[int] = []
    try:
        with gzip.open(expression_path, "rt", encoding="utf-8-sig", newline="") as source:
            reader = csv.reader(source, delimiter="\t")
            header = next(reader, None)
            if header is None or len(header) < 2:
                raise ValueError("Expression source must have one gene column and sample columns")
            normalized_header = [value.strip() for value in header]
            if any(
                original != normalized
                for original, normalized in zip(header, normalized_header, strict=True)
            ):
                raise ValueError("Expression source contains whitespace-padded column identifiers")
            if normalized_header[0].casefold() != "sample":
                raise ValueError(
                    f"Expression source first column must be 'sample', found {header[0]!r}"
                )
            if len(normalized_header) != len(set(normalized_header)):
                raise ValueError("Expression source contains duplicate column identifiers")
            sample_header = normalized_header[1:]
            source_sample_count = len(sample_header)
            index_by_sample = {
                sample_id: index + 1 for index, sample_id in enumerate(sample_header)
            }
            missing_samples = sorted(set(sample_ids) - set(index_by_sample))
            if missing_samples:
                raise ValueError(
                    f"Cohort samples missing from expression source: {missing_samples[:10]}"
                )
            sample_indices = [index_by_sample[sample_id] for sample_id in sample_ids]
            noncohort_source_sample_count = source_sample_count - len(sample_ids)
            if noncohort_source_sample_count < 0:
                raise ValueError("Cohort sample count exceeds expression source sample count")

            with temp_matrix_path.open("wb") as raw_output:
                with gzip.GzipFile(
                    filename="",
                    fileobj=raw_output,
                    mode="wb",
                    compresslevel=6,
                    mtime=0,
                ) as compressed_output:
                    with io.TextIOWrapper(
                        compressed_output,
                        encoding="utf-8",
                        newline="",
                    ) as text_output:
                        writer = csv.writer(
                            text_output,
                            delimiter="\t",
                            lineterminator="\n",
                        )
                        writer.writerow(["gene_id", *sample_ids])
                        for row_number, row in enumerate(reader, start=2):
                            if len(row) != len(normalized_header):
                                raise ValueError(
                                    f"Expression row {row_number} has {len(row)} fields; "
                                    f"expected {len(normalized_header)}"
                                )
                            gene_id = row[0]
                            if not gene_id:
                                raise ValueError(f"Expression row {row_number} has a blank gene ID")
                            if gene_id != gene_id.strip():
                                raise ValueError(
                                    f"Expression row {row_number} has a whitespace-padded gene ID"
                                )
                            if gene_id in seen_source_genes:
                                raise ValueError(f"Duplicate source gene ID: {gene_id}")
                            seen_source_genes.add(gene_id)
                            source_gene_count += 1
                            if (
                                expected_source_gene_count is not None
                                and source_gene_count > expected_source_gene_count
                            ):
                                raise ValueError("Expression source has more genes than configured")
                            if source_gene_count > len(all_gene_rows):
                                raise ValueError(
                                    "Expression source has more rows than Phase 1D summary"
                                )
                            gene_summary_row = all_gene_rows[source_gene_count - 1]
                            if gene_id != gene_summary_row["gene_id"]:
                                raise ValueError(
                                    "Expression source gene order differs from the Phase 1D "
                                    f"gene summary at row {source_gene_count}: {gene_id!r} != "
                                    f"{gene_summary_row['gene_id']!r}"
                                )
                            gene_summary_id = gene_summary_row["gene_id"]
                            n_above = 0
                            selected_values: list[str] = []
                            for sample_id, column_index in zip(
                                sample_ids, sample_indices, strict=True
                            ):
                                source_token = row[column_index]
                                numeric_token = source_token.strip()
                                if numeric_token.casefold() in MISSING_TOKENS:
                                    raise ValueError(
                                        "Missing expression value for selected sample "
                                        f"{sample_id}, gene {gene_id}"
                                    )
                                try:
                                    value = float(numeric_token)
                                except ValueError as exc:
                                    raise ValueError(
                                        f"Non-numeric expression for sample {sample_id}, "
                                        f"gene {gene_id}: {source_token!r}"
                                    ) from exc
                                if not math.isfinite(value):
                                    raise ValueError(
                                        f"Non-finite expression for sample {sample_id}, "
                                        f"gene {gene_id}: {source_token!r}"
                                    )
                                n_above += value > expression_threshold
                                selected_values.append(source_token)
                            recorded_count = _as_int(
                                gene_summary_row["n_samples_above_selected_threshold"],
                                field="n_samples_above_selected_threshold",
                                row_id=gene_summary_id,
                            )
                            if n_above != recorded_count:
                                raise ValueError(
                                    f"Source prevalence count for {gene_id} is {n_above}, "
                                    f"but Phase 1D records {recorded_count}"
                                )
                            retain = n_above >= required_samples
                            recorded_flag = _as_bool(
                                gene_summary_row["retained_by_prevalence_filter"],
                                field="retained_by_prevalence_filter",
                                row_id=gene_id,
                            )
                            if retain != recorded_flag:
                                raise ValueError(
                                    "Source prevalence filter disagrees with Phase 1D "
                                    f"for {gene_id}"
                                )
                            if retain:
                                if selected_index >= len(selected_gene_ids):
                                    raise ValueError(
                                        "Source contains more passing genes than summary"
                                    )
                                if gene_id != selected_gene_ids[selected_index]:
                                    raise ValueError(
                                        "Retained gene order differs from the Phase 1D "
                                        "selected gene list"
                                    )
                                writer.writerow([gene_id, *selected_values])
                                selected_index += 1

        if (
            expected_source_gene_count is not None
            and source_gene_count != expected_source_gene_count
        ):
            raise ValueError(
                f"Expression source has {source_gene_count} genes; "
                f"expected {expected_source_gene_count}"
            )
        if source_gene_count != len(all_gene_rows):
            raise ValueError(
                f"Expression source has {source_gene_count} genes; "
                f"Phase 1D summary has {len(all_gene_rows)}"
            )
        if selected_index != len(selected_gene_ids):
            raise ValueError(
                f"Wrote {selected_index} retained genes; expected {len(selected_gene_ids)}"
            )
        temp_matrix_path.replace(output_path)
    finally:
        if temp_matrix_path.exists():
            temp_matrix_path.unlink()

    _verify_output_dimensions(
        output_path,
        sample_ids=sample_ids,
        gene_ids=selected_gene_ids,
    )
    _write_tsv(gene_list_path, selected_gene_rows, GENE_LIST_FIELDS)
    matrix_checksum = sha256_file(output_path)
    gene_list_checksum = sha256_file(gene_list_path)

    input_paths = [
        expression_path,
        cohort_path,
        gene_summary_path,
        sample_metrics_path,
        Path(config_path),
    ]
    inputs = [_input_record(path.resolve(), project_root) for path in input_paths]
    manifest = {
        "manifest_version": 1,
        "created_date": datetime.now(ZoneInfo("Asia/Singapore")).date().isoformat(),
        "inputs": inputs,
        "expression_units": expression_scale,
        "source_matrix": {
            "path": _manifest_path(expression_path, project_root),
            "gene_rows": source_gene_count,
            "sample_columns": source_sample_count,
            "cohort_sample_columns_used": len(sample_ids),
            "noncohort_sample_columns_omitted": noncohort_source_sample_count,
        },
        "matrix": {
            "path": output_path.relative_to(project_root).as_posix(),
            "sha256": matrix_checksum,
            "format": "UTF-8 TSV compressed with gzip",
            "orientation": "genes are rows; samples are columns",
            "first_column": "gene_id",
            "sample_order_source": "classification_cohort.tsv row order",
            "gene_order_source": (
                "source expression row order, validated against qc_gene_summary.tsv"
            ),
            "additional_transformation": "none; selected source numeric tokens copied verbatim",
            "shape": {
                "gene_rows": len(selected_gene_ids),
                "sample_columns": len(sample_ids),
                "total_columns_including_gene_id": len(sample_ids) + 1,
            },
        },
        "gene_list": {
            "path": gene_list_path.relative_to(project_root).as_posix(),
            "sha256": gene_list_checksum,
            "gene_order_is_source_order": True,
        },
        "selected_filter": {
            "rule": (
                f"finite expression > {expression_threshold:g} in at least "
                f"{minimum_sample_fraction:.0%} of cohort samples"
            ),
            "expression_threshold_strictly_greater_than": expression_threshold,
            "minimum_sample_fraction": minimum_sample_fraction,
            "required_sample_count": required_samples,
            "sample_count": len(sample_ids),
            "source_gene_count": source_gene_count,
            "retained_gene_count": len(selected_gene_ids),
        },
        "configuration": dict(config_snapshot),
        "sample_order": sample_ids,
        "gene_order": selected_gene_ids,
        "tukey_outlier_samples_flagged_but_retained": flagged_outliers,
        "clinical_or_pam50_columns_in_expression_matrix": False,
        "splits_or_modeling_performed": False,
    }
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    temp_manifest = manifest_path.with_name(manifest_path.name + ".tmp")
    try:
        with temp_manifest.open("w", encoding="utf-8", newline="\n") as handle:
            json.dump(manifest, handle, indent=2, sort_keys=True)
            handle.write("\n")
        temp_manifest.replace(manifest_path)
    finally:
        if temp_manifest.exists():
            temp_manifest.unlink()

    return MLMatrixBuildResult(
        gene_ids=tuple(selected_gene_ids),
        sample_ids=tuple(sample_ids),
        source_gene_count=source_gene_count,
        source_sample_count=source_sample_count,
        noncohort_source_sample_count=noncohort_source_sample_count,
        required_sample_count=required_samples,
        flagged_outlier_sample_ids=tuple(flagged_outliers),
        matrix_path=output_path,
        gene_list_path=gene_list_path,
        manifest_path=manifest_path,
        matrix_sha256=matrix_checksum,
    )
