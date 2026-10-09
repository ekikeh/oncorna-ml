"""Build the local-only all-gene view and fit a training-only prevalence filter.

This command reads split_v1 and Phase 1 artifacts but never regenerates them.
It creates one new versioned directory under Git-ignored data/processed/.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
import tempfile
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import sklearn
import yaml

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_PATH = PROJECT_ROOT / "src"
if str(SRC_PATH) not in sys.path:
    sys.path.insert(0, str(SRC_PATH))

from oncorna.preprocessing import (  # noqa: E402
    GenePrevalenceFilter,
    identity_digest,
    partition_frames,
    read_all_gene_view,
    sha256_file,
    validate_frozen_partitions,
    verify_split_provenance,
)
from oncorna.qc import read_classification_cohort  # noqa: E402

DEFAULT_CONFIG = PROJECT_ROOT / "configs" / "preprocessing_v1.yaml"


def _project_path(value: str) -> Path:
    path = Path(value)
    return path if path.is_absolute() else PROJECT_ROOT / path


def _write_json(path: Path, value: Any) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(value, handle, ensure_ascii=False, indent=2, sort_keys=True)
        handle.write("\n")


def _sequence_digest(values: tuple[str, ...]) -> str:
    payload = json.dumps(values, ensure_ascii=False, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _check_source_manifest(source_path: Path, expected_sha256: str) -> None:
    manifest_path = PROJECT_ROOT / "data" / "MANIFEST.sha256"
    relative_source = source_path.relative_to(PROJECT_ROOT).as_posix()
    entries = [line.split() for line in manifest_path.read_text(encoding="utf-8").splitlines()]
    matches = [fields[0] for fields in entries if len(fields) == 2 and fields[1] == relative_source]
    if matches != [expected_sha256]:
        raise ValueError("Tracked source checksum manifest disagrees with preprocessing config")
    if sha256_file(source_path) != expected_sha256:
        raise ValueError("Raw expression source SHA-256 differs from its recorded checksum")


def prepare(config_path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    """Validate frozen inputs, fit on train only, and publish new local artifacts."""
    config_path = Path(config_path).resolve()
    with config_path.open("r", encoding="utf-8") as handle:
        config = yaml.safe_load(handle)
    if config["version"] != "preprocessing_v1":
        raise ValueError("Expected preprocessing_v1 configuration")

    expression_path = _project_path(config["expression_path"]).resolve()
    cohort_path = _project_path(config["cohort_path"]).resolve()
    split_path = _project_path(config["split_path"]).resolve()
    historical_config_path = _project_path(config["historical_config_path"]).resolve()
    historical_matrix_path = _project_path(config["historical_matrix_path"]).resolve()
    historical_matrix_manifest_path = _project_path(
        config["historical_matrix_manifest_path"]
    ).resolve()
    output_dir = _project_path(config["output_dir"]).resolve()
    processed_dir = (PROJECT_ROOT / "data" / "processed").resolve()
    if not output_dir.is_relative_to(processed_dir) or output_dir == processed_dir:
        raise ValueError("New preprocessing output must be under Git-ignored data/processed/")
    if output_dir.exists():
        raise FileExistsError(f"Refusing to overwrite existing preprocessing output: {output_dir}")

    expected_source_sha256 = config["expected_source_sha256"]
    _check_source_manifest(expression_path, expected_source_sha256)
    split_sha256 = sha256_file(split_path)
    if split_sha256 != config["expected_split_sha256"]:
        raise ValueError("Frozen split SHA-256 differs from the approved checkpoint")
    with split_path.open("r", encoding="utf-8") as handle:
        split = json.load(handle)
    cohort_records = read_classification_cohort(
        cohort_path, expected_sample_count=int(config["expected_cohort_size"])
    )
    partitions = validate_frozen_partitions(cohort_records, split)
    provenance = verify_split_provenance(
        split,
        cohort_path=cohort_path,
        historical_config_path=historical_config_path,
        historical_matrix_path=historical_matrix_path,
        historical_matrix_manifest_path=historical_matrix_manifest_path,
    )
    provenance["historical_config"]["path"] = historical_config_path.relative_to(
        PROJECT_ROOT
    ).as_posix()

    cohort_sample_ids = tuple(row["sample_id"] for row in cohort_records)
    view = read_all_gene_view(
        expression_path,
        cohort_sample_ids,
        expected_gene_count=int(config["expected_gene_count"]),
    )
    if view.sample_ids != cohort_sample_ids:
        raise ValueError("All-gene sample order differs from the approved cohort")
    missing_values = int(np.isnan(view.values).sum())
    infinite_values = int(np.isinf(view.values).sum())
    if infinite_values:
        raise ValueError("All-gene view contains infinite expression")

    frames = partition_frames(view, partitions)
    training_pairs = partitions["train"]
    if len(training_pairs) != int(config["expected_training_size"]):
        raise ValueError("Frozen training partition size differs from approved count")
    gene_filter = GenePrevalenceFilter(
        expression_gt=float(config["expression_gt"]),
        minimum_sample_fraction=float(config["minimum_sample_fraction"]),
    ).fit(frames["train"])
    if gene_filter.fitting_sample_ids_ != tuple(sample for _, sample in training_pairs):
        raise ValueError("Filter fit sample IDs differ from frozen training assignments")
    if gene_filter.required_count_ != 102:
        raise ValueError("Approved full-training prevalence count must be 102")
    if not gene_filter.selected_gene_ids_:
        raise ValueError("Training-only filter retained no genes")
    for name in ("train", "validation", "test"):
        transformed = gene_filter.transform(frames[name])
        if tuple(transformed.columns) != gene_filter.selected_gene_ids_:
            raise ValueError(f"{name} transformed gene order differs from training")
        if tuple(transformed.index) != tuple(sample for _, sample in partitions[name]):
            raise ValueError(f"{name} transformed sample order differs from the split")

    processed_dir.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".preprocessing-v1-", dir=processed_dir) as temp:
        temp_dir = Path(temp)
        matrix_path = temp_dir / "all_gene_expression.npy"
        schema_path = temp_dir / "all_gene_schema.json"
        gene_list_path = temp_dir / "training_gene_list.tsv"
        manifest_path = temp_dir / "manifest.json"
        with matrix_path.open("wb") as handle:
            np.save(handle, view.values, allow_pickle=False)
        _write_json(
            schema_path,
            {
                "version": "preprocessing_v1",
                "orientation": "samples_by_genes",
                "expression_units": "log2(normalized_count + 1)",
                "dtype": "float64",
                "shape": list(view.values.shape),
                "sample_ids": list(view.sample_ids),
                "gene_ids": list(view.gene_ids),
            },
        )
        with gene_list_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(("gene_id", "training_samples_gt_1"))
            for gene_id, count, keep in zip(
                view.gene_ids,
                gene_filter.prevalence_counts_,
                gene_filter.support_mask_,
                strict=True,
            ):
                if keep:
                    writer.writerow((gene_id, int(count)))
        manifest = {
            "version": "preprocessing_v1",
            "source": {
                "path": expression_path.relative_to(PROJECT_ROOT).as_posix(),
                "sha256": expected_source_sha256,
                "source_sample_count": view.source_sample_count,
                "gene_count": len(view.gene_ids),
            },
            "cohort": {
                "path": cohort_path.relative_to(PROJECT_ROOT).as_posix(),
                "sha256": provenance["cohort_sha256"],
                "sample_count": len(view.sample_ids),
                "sample_order_sha256": _sequence_digest(view.sample_ids),
            },
            "split": {
                "path": split_path.relative_to(PROJECT_ROOT).as_posix(),
                "sha256": split_sha256,
                "partition_counts": {name: len(pairs) for name, pairs in partitions.items()},
                "partition_identity_sha256": {
                    name: identity_digest(pairs) for name, pairs in partitions.items()
                },
            },
            "historical_provenance": provenance,
            "preprocessing_config": {
                "path": config_path.relative_to(PROJECT_ROOT).as_posix(),
                "sha256": sha256_file(config_path),
            },
            "all_gene_view": {
                "path": f"{config['output_dir']}/all_gene_expression.npy",
                "sha256": sha256_file(matrix_path),
                "schema_path": f"{config['output_dir']}/all_gene_schema.json",
                "schema_sha256": sha256_file(schema_path),
                "shape": list(view.values.shape),
                "gene_order_sha256": _sequence_digest(view.gene_ids),
                "missing_values": missing_values,
                "infinite_values": infinite_values,
            },
            "filter": {
                "fit_partition": "train",
                "fit_patient_sample_sha256": identity_digest(training_pairs),
                "fit_patient_count": len(training_pairs),
                "fit_sample_count": gene_filter.n_samples_fit_,
                "expression_gt": gene_filter.expression_gt,
                "minimum_sample_fraction": gene_filter.minimum_sample_fraction,
                "required_sample_count": gene_filter.required_count_,
                "missing_values_in_fit": gene_filter.missing_value_count_,
                "infinite_values_in_fit": gene_filter.infinite_value_count_,
                "retained_gene_count": len(gene_filter.selected_gene_ids_),
                "gene_list_path": f"{config['output_dir']}/training_gene_list.tsv",
                "gene_list_sha256": sha256_file(gene_list_path),
            },
            "transform_check": {
                "partitions": ["train", "validation", "test"],
                "reused_training_gene_mask": True,
                "output_gene_order_sha256": _sequence_digest(gene_filter.selected_gene_ids_),
            },
            "software": {
                "python": sys.version.split()[0],
                "numpy": np.__version__,
                "pandas": pd.__version__,
                "scikit_learn": sklearn.__version__,
            },
        }
        _write_json(manifest_path, manifest)
        if output_dir.exists():
            raise FileExistsError(
                f"Refusing to overwrite existing preprocessing output: {output_dir}"
            )
        temp_dir.rename(output_dir)
    if sha256_file(split_path) != split_sha256:
        raise ValueError("Frozen split changed during preprocessing")
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    args = parser.parse_args()
    manifest = prepare(args.config)
    print(f"All-gene view: {manifest['source']['gene_count']:,} genes x 844 samples")
    print(f"Missing values: {manifest['all_gene_view']['missing_values']:,}")
    print(f"Infinite values: {manifest['all_gene_view']['infinite_values']:,}")
    print(
        "Training-only filter: "
        f"{manifest['filter']['fit_sample_count']} samples, "
        f"> {manifest['filter']['expression_gt']:g} in "
        f">= {manifest['filter']['required_sample_count']} samples"
    )
    print(f"Retained genes: {manifest['filter']['retained_gene_count']:,}")
    print(f"Gene-list SHA-256: {manifest['filter']['gene_list_sha256']}")
    print(f"Manifest: {manifest['filter']['gene_list_path'].rsplit('/', 1)[0]}/manifest.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
