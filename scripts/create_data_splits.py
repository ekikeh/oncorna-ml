"""Create the frozen Phase 1F patient-level train/validation/test split.

Run from the repository root with:

    python scripts/create_data_splits.py

The script uses cohort IDs and normalized PAM50 labels only. It checks the
expression matrix header but never reads expression values.
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

from oncorna.splits import generate_split_manifest, write_frozen_split_manifest  # noqa: E402

DEFAULT_CONFIG = PROJECT_ROOT / "configs" / "default.yaml"


def _load_config(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        config = yaml.safe_load(handle)
    if not isinstance(config, dict):
        raise ValueError(f"Configuration file must contain a YAML mapping: {path}")
    for section in ("qc", "ml_matrix", "splitting"):
        if not isinstance(config.get(section), dict):
            raise ValueError(f"Configuration file has no {section} section: {path}")
    return config


def _project_path(value: str) -> Path:
    path = Path(value)
    return path if path.is_absolute() else PROJECT_ROOT / path


def create_splits(config_path: Path = DEFAULT_CONFIG) -> tuple[dict[str, Any], bool]:
    """Generate split metadata and create it once, without replacing split_v1."""
    config_path = Path(config_path)
    if not config_path.is_absolute():
        config_path = PROJECT_ROOT / config_path
    config_path = config_path.resolve()
    config = _load_config(config_path)
    qc = config["qc"]
    matrix = config["ml_matrix"]
    splitting = config["splitting"]

    if splitting.get("split_unit") != "patient_id":
        raise ValueError("Phase 1F split_unit must be patient_id")
    if splitting.get("stratification_label") != "pam50_normalized_label":
        raise ValueError("Phase 1F stratification_label must be pam50_normalized_label")
    output_path = _project_path(splitting["output_path"]).resolve()
    processed_dir = (PROJECT_ROOT / "data" / "processed").resolve()
    expected_output_path = (processed_dir / "split_v1.json").resolve()
    if output_path != expected_output_path:
        raise ValueError("Phase 1F output_path must be data/processed/split_v1.json")
    if not output_path.is_relative_to(processed_dir):
        raise ValueError("Phase 1F output must remain under Git-ignored data/processed/")

    manifest = generate_split_manifest(
        cohort_path=_project_path(qc["cohort_path"]),
        expression_matrix_path=_project_path(matrix["output_path"]),
        matrix_manifest_path=_project_path(matrix["manifest_path"]),
        config_path=config_path,
        project_root=PROJECT_ROOT,
        expected_patient_count=int(qc["expected_cohort_samples"]),
        expected_gene_count=int(matrix["expected_gene_count"]),
        expression_scale=str(qc["expression_scale"]),
        random_seed=int(config["random_seed"]),
        proportions={key: float(value) for key, value in splitting["proportions"].items()},
        rounding_tie_order=list(splitting["rounding_tie_order"]),
        test_random_state_offset=int(splitting["test_random_state_offset"]),
        validation_random_state_offset=int(splitting["validation_random_state_offset"]),
    )
    manifest["output_path"] = output_path.relative_to(PROJECT_ROOT).as_posix()
    created = write_frozen_split_manifest(manifest, output_path)
    return manifest, created


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config",
        type=Path,
        default=DEFAULT_CONFIG,
        help=f"configuration YAML (default: {DEFAULT_CONFIG.relative_to(PROJECT_ROOT)})",
    )
    args = parser.parse_args()
    config_path = Path(args.config)
    if not config_path.is_absolute():
        config_path = PROJECT_ROOT / config_path
    config_path = config_path.resolve()
    manifest, created = create_splits(config_path)
    output_path = _project_path(_load_config(config_path)["splitting"]["output_path"])
    action = "Created and froze" if created else "Verified existing frozen"
    print(f"{action} {manifest['split_id']}")
    print(f"Source cohort SHA-256: {manifest['source_cohort']['sha256']}")
    print(f"Random seed: {manifest['random_seed']}")
    for partition in ("train", "validation", "test"):
        part = manifest["partitions"][partition]
        classes = ", ".join(f"{label}={count}" for label, count in part["pam50_counts"].items())
        print(
            f"{partition}: {part['n_patients']} patients / {part['n_samples']} samples; {classes}"
        )
    print(f"Validation checks passed: {manifest['validation_results']['passed']}")
    print(f"Output: {output_path.relative_to(PROJECT_ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
