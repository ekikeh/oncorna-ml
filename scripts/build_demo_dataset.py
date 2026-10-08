"""Build the small Phase 1G demo tables from local frozen Phase 1C/1E/1F inputs.

Run from the repository root with:

    python scripts/build_demo_dataset.py

The command does not create or modify the frozen split. It only reads the
approved cohort, split_v1, processed matrix, and retained-gene list, then writes
small outputs under data/demo/.
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

from oncorna.demo import build_demo_dataset  # noqa: E402

DEFAULT_CONFIG = PROJECT_ROOT / "configs" / "default.yaml"
DEFAULT_DEMO_CONFIG = PROJECT_ROOT / "configs" / "demo.yaml"


def _load_yaml(path: Path) -> dict[str, Any]:
    with Path(path).open("r", encoding="utf-8") as handle:
        config = yaml.safe_load(handle)
    if not isinstance(config, dict):
        raise ValueError(f"Configuration file must contain a YAML mapping: {path}")
    return config


def _project_path(value: str) -> Path:
    path = Path(value)
    return path if path.is_absolute() else PROJECT_ROOT / path


def create_demo(
    config_path: Path = DEFAULT_CONFIG,
    demo_config_path: Path = DEFAULT_DEMO_CONFIG,
):
    """Build the demo using the existing project configuration and Phase 1G settings."""
    config_path = Path(config_path)
    if not config_path.is_absolute():
        config_path = PROJECT_ROOT / config_path
    config_path = config_path.resolve()
    demo_config_path = Path(demo_config_path)
    if not demo_config_path.is_absolute():
        demo_config_path = PROJECT_ROOT / demo_config_path
    demo_config_path = demo_config_path.resolve()

    config = _load_yaml(config_path)
    demo = _load_yaml(demo_config_path)
    for section in ("qc", "ml_matrix", "splitting"):
        if not isinstance(config.get(section), dict):
            raise ValueError(f"Project configuration has no {section} section: {config_path}")
    for field in (
        "partition",
        "samples_per_class",
        "gene_count",
        "output_expression_path",
        "output_metadata_path",
        "output_manifest_path",
    ):
        if field not in demo:
            raise ValueError(f"Demo configuration is missing {field!r}: {demo_config_path}")
    if demo["partition"] != "train":
        raise ValueError("Phase 1G demo samples must come only from the frozen training partition")

    qc = config["qc"]
    ml_matrix = config["ml_matrix"]
    splitting = config["splitting"]
    return build_demo_dataset(
        cohort_path=_project_path(qc["cohort_path"]),
        split_manifest_path=_project_path(splitting["output_path"]),
        expression_matrix_path=_project_path(ml_matrix["output_path"]),
        gene_list_path=_project_path(ml_matrix["gene_list_path"]),
        matrix_manifest_path=_project_path(ml_matrix["manifest_path"]),
        project_config_path=config_path,
        demo_config_path=demo_config_path,
        project_root=PROJECT_ROOT,
        expected_sample_count=int(qc["expected_cohort_samples"]),
        expected_gene_count=int(ml_matrix["expected_gene_count"]),
        expression_scale=str(qc["expression_scale"]),
        random_seed=int(config["random_seed"]),
        samples_per_class=int(demo["samples_per_class"]),
        demo_gene_count=int(demo["gene_count"]),
        output_expression_path=_project_path(demo["output_expression_path"]),
        output_metadata_path=_project_path(demo["output_metadata_path"]),
        output_manifest_path=_project_path(demo["output_manifest_path"]),
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config",
        type=Path,
        default=DEFAULT_CONFIG,
        help=f"project configuration (default: {DEFAULT_CONFIG.relative_to(PROJECT_ROOT)})",
    )
    parser.add_argument(
        "--demo-config",
        type=Path,
        default=DEFAULT_DEMO_CONFIG,
        help=f"Phase 1G settings (default: {DEFAULT_DEMO_CONFIG.relative_to(PROJECT_ROOT)})",
    )
    args = parser.parse_args()
    result = create_demo(args.config, args.demo_config)
    print("Built Phase 1G demo dataset")
    print(f"Samples: {len(result.sample_ids)}; genes: {len(result.gene_ids)}")
    print(f"PAM50 counts: {result.class_counts}")
    print(f"Expression: {result.expression_path.relative_to(PROJECT_ROOT)}")
    print(f"Metadata: {result.metadata_path.relative_to(PROJECT_ROOT)}")
    print(f"Manifest: {result.manifest_path.relative_to(PROJECT_ROOT)}")
    print("Frozen split_v1: read-only and checksum-verified")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
