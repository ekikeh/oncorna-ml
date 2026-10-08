"""Build the Phase 1E expression-only matrix from frozen Phase 1C/1D artifacts.

Run from the repository root with:

    python scripts/build_ml_matrix.py

The full matrix, selected-gene list, and manifest are written under
``data/processed/``, which is Git-ignored.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Any

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_PATH = PROJECT_ROOT / "src"
if str(SRC_PATH) not in sys.path:
    sys.path.insert(0, str(SRC_PATH))

import yaml  # noqa: E402

from oncorna.ml_matrix import MLMatrixBuildResult, build_ml_matrix  # noqa: E402

DEFAULT_CONFIG = PROJECT_ROOT / "configs" / "default.yaml"


def _load_config(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        config = yaml.safe_load(handle)
    if not isinstance(config, dict):
        raise ValueError(f"Configuration file must contain a YAML mapping: {path}")
    for section in ("qc", "ml_matrix"):
        if not isinstance(config.get(section), dict):
            raise ValueError(f"Configuration file has no {section} section: {path}")
    return config


def _project_path(root: Path, value: str) -> Path:
    path = Path(value)
    return path if path.is_absolute() else root / path


def _manifest_configuration(config: dict[str, Any]) -> dict[str, Any]:
    qc = config["qc"]
    matrix = config["ml_matrix"]
    return {
        "qc": {
            "expected_cohort_samples": qc["expected_cohort_samples"],
            "expression_scale": qc["expression_scale"],
            "gene_filter": dict(qc["gene_filter"]),
            "sample_filter": dict(qc["sample_filter"]),
            "outlier_detection": dict(qc["outlier_detection"]),
        },
        "ml_matrix": dict(matrix),
    }


def build_from_config(config_path: Path = DEFAULT_CONFIG) -> MLMatrixBuildResult:
    """Load configured Phase 1C/1D inputs and build the local-only matrix."""
    config_path = Path(config_path)
    if not config_path.is_absolute():
        config_path = PROJECT_ROOT / config_path
    config_path = config_path.resolve()
    config = _load_config(config_path)
    qc = config["qc"]
    matrix = config["ml_matrix"]

    required_matrix_keys = (
        "expected_sample_count",
        "expected_source_gene_count",
        "expected_gene_count",
        "output_path",
        "gene_list_path",
        "manifest_path",
        "expression_transform",
    )
    missing = [key for key in required_matrix_keys if key not in matrix]
    if missing:
        raise ValueError(f"ml_matrix config is missing required settings: {missing}")

    if int(qc["expected_cohort_samples"]) != int(matrix["expected_sample_count"]):
        raise ValueError(
            "qc.expected_cohort_samples and ml_matrix.expected_sample_count must agree"
        )

    cohort_path = _project_path(PROJECT_ROOT, qc["cohort_path"])
    expression_path = _project_path(PROJECT_ROOT, qc["expression_path"])
    gene_summary_path = _project_path(PROJECT_ROOT, qc["gene_summary_path"])
    sample_metrics_path = _project_path(PROJECT_ROOT, qc["sample_metrics_path"])
    return build_ml_matrix(
        expression_path=expression_path,
        cohort_path=cohort_path,
        gene_summary_path=gene_summary_path,
        sample_metrics_path=sample_metrics_path,
        config_path=config_path,
        output_path=_project_path(PROJECT_ROOT, matrix["output_path"]),
        gene_list_path=_project_path(PROJECT_ROOT, matrix["gene_list_path"]),
        manifest_path=_project_path(PROJECT_ROOT, matrix["manifest_path"]),
        project_root=PROJECT_ROOT,
        expression_scale=qc["expression_scale"],
        expression_threshold=float(qc["gene_filter"]["expression_gt"]),
        minimum_sample_fraction=float(qc["gene_filter"]["minimum_sample_fraction"]),
        expected_sample_count=int(matrix["expected_sample_count"]),
        expected_gene_count=int(matrix["expected_gene_count"]),
        expected_source_gene_count=int(matrix["expected_source_gene_count"]),
        expression_transform=str(matrix["expression_transform"]),
        config_snapshot=_manifest_configuration(config),
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config",
        type=Path,
        default=DEFAULT_CONFIG,
        help=f"configuration YAML (default: {DEFAULT_CONFIG.relative_to(PROJECT_ROOT)})",
    )
    args = parser.parse_args()
    result = build_from_config(args.config)
    print(
        "Built Phase 1E matrix: "
        f"{len(result.gene_ids):,} genes x {len(result.sample_ids):,} samples"
    )
    print(
        f"Required prevalence: > configured threshold in >= {result.required_sample_count} samples"
    )
    print(
        "Source dimensions: "
        f"{result.source_gene_count:,} genes x "
        f"{result.source_sample_count:,} samples"
    )
    print(f"Non-cohort source samples omitted: {result.noncohort_source_sample_count:,}")
    print(f"Flagged Tukey outliers retained: {len(result.flagged_outlier_sample_ids)}")
    print(f"Matrix: {result.matrix_path.relative_to(PROJECT_ROOT)}")
    print(f"Gene list: {result.gene_list_path.relative_to(PROJECT_ROOT)}")
    print(f"Manifest: {result.manifest_path.relative_to(PROJECT_ROOT)}")
    print(f"Matrix SHA-256: {result.matrix_sha256}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
