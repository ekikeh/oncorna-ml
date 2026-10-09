"""Fit the approved 506-patient model and evaluate frozen validation exactly once."""

from __future__ import annotations

import argparse
import csv
import json
import platform
import subprocess
import sys
import traceback
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import joblib
import numpy as np
import pandas as pd
import sklearn
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from oncorna.preprocessing import sha256_file  # noqa: E402
from oncorna.validation import (  # noqa: E402
    fit_and_evaluate,
    load_validation_inputs,
    verify_phase2b,
)


def write_json(path: Path, value: Any) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(value, handle, indent=2, sort_keys=True, ensure_ascii=False)
        handle.write("\n")


def git_value(*args: str) -> str:
    return subprocess.check_output(["git", *args], cwd=ROOT, text=True).strip()


def preflight(config_path: Path) -> tuple[dict[str, Any], dict[str, Any], Any, Path]:
    approval = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    if (
        approval["version"] != "validation_v1"
        or approval["approved_commit"] != "f58abd350d0d49295249c29f6fbf898620987481"
        or approval["selected_C"] != 1.0
        or (approval["training_size"], approval["validation_size"], approval["test_size"])
        != (506, 169, 169)
    ):
        raise ValueError("Phase 2C approval configuration differs from authorized protocol")
    if git_value("rev-parse", "HEAD") != approval["approved_commit"]:
        raise ValueError("Git HEAD differs from approved Phase 2B commit")
    if git_value("branch", "--show-current") != "main":
        raise ValueError("Phase 2C requires the approved main checkout")
    output = (ROOT / approval["output_dir"]).resolve()
    processed = (ROOT / "data/processed").resolve()
    if not output.is_relative_to(processed) or output == processed or output.exists():
        raise ValueError("Validation output must be a new local-only versioned directory")
    modeling_path = ROOT / approval["modeling_config"]
    if sha256_file(modeling_path) != approval["modeling_config_sha256"]:
        raise ValueError("Approved Phase 2B configuration checksum mismatch")
    modeling = yaml.safe_load(modeling_path.read_text(encoding="utf-8"))
    if (
        modeling["labels"]
        != ["Luminal A", "Luminal B", "Basal-like", "HER2-enriched", "Normal-like"]
        or modeling["expression_gt"] != 1.0
        or modeling["minimum_sample_fraction"] != 0.20
        or modeling["solver"] != "lbfgs"
        or modeling["penalty"] != "l2"
        or modeling["class_weight"] is not None
        or modeling["fit_intercept"] is not True
        or modeling["max_iter"] != 1000
        or modeling["tol"] != 0.0001
    ):
        raise ValueError("Modeling configuration differs from approved Phase 2C model")
    phase2b_hash = verify_phase2b(ROOT, approval, modeling)
    inputs = load_validation_inputs(ROOT, approval, modeling, phase2b_hash)
    return approval, modeling, inputs, output


def run(config_path: Path) -> dict[str, Any]:
    config_path = config_path.resolve()
    approval, modeling, inputs, output = preflight(config_path)
    dirty_at_run = git_value("status", "--porcelain")
    started_at = datetime.now(timezone.utc).isoformat()
    output.mkdir(parents=True)
    try:
        pipeline, results, rows = fit_and_evaluate(inputs, modeling, approval["selected_C"])
        joblib.dump(pipeline, output / "pipeline.joblib")
        write_json(output / "training_manifest.json", results["training"])
        write_json(output / "validation_metrics.json", results["validation"])
        write_json(
            output / "confusion_matrices.json",
            {
                "class_order": modeling["labels"],
                "raw": results["validation"]["confusion_matrix"],
                "row_normalized": results["validation"]["row_normalized_confusion_matrix"],
            },
        )
        with (output / "validation_predictions.tsv").open(
            "w", encoding="utf-8", newline=""
        ) as handle:
            fields = [
                "patient_id",
                "sample_id",
                "true_label",
                "predicted_label",
                *(f"probability_{label}" for label in modeling["labels"]),
            ]
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        gene_ids = pipeline.named_steps["gene_filter"].selected_gene_ids_
        with (output / "training_gene_list.tsv").open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["gene_id"])
            writer.writerows((gene,) for gene in gene_ids)
        output_names = (
            "pipeline.joblib",
            "training_manifest.json",
            "validation_metrics.json",
            "confusion_matrices.json",
            "validation_predictions.tsv",
            "training_gene_list.tsv",
        )
        manifest = {
            "version": "validation_v1",
            "scope": "single_frozen_validation_evaluation",
            "evaluation_count": 1,
            "started_at_utc": started_at,
            "finished_at_utc": datetime.now(timezone.utc).isoformat(),
            "approved_git_commit": approval["approved_commit"],
            "git_commit_at_run": git_value("rev-parse", "HEAD"),
            "git_dirty_at_run": bool(dirty_at_run),
            "config_sha256": sha256_file(config_path),
            "modeling_config_sha256": approval["modeling_config_sha256"],
            "phase2b_manifest_sha256": inputs.phase2b_manifest_sha256,
            "input_hashes": inputs.training.input_hashes,
            "partition_counts": {
                "train": len(inputs.training.y),
                "validation": len(inputs.y_validation),
                "test_unaccessed": approval["test_size"],
            },
            "class_order": modeling["labels"],
            "selected_C": approval["selected_C"],
            "software": {
                "python": platform.python_version(),
                "numpy": np.__version__,
                "pandas": pd.__version__,
                "scikit_learn": sklearn.__version__,
                "joblib": joblib.__version__,
            },
            "output_hashes": {name: sha256_file(output / name) for name in output_names},
        }
        write_json(output / "manifest.json", manifest)
    except Exception as exc:
        write_json(
            output / "failure.json",
            {
                "error_type": type(exc).__name__,
                "message": str(exc),
                "traceback": traceback.format_exc(),
                "started_at_utc": started_at,
            },
        )
        raise
    print(f"Training genes retained: {results['training']['retained_gene_count']}")
    print(f"Validation macro F1: {results['validation']['macro_f1']:.4f}")
    print(f"Manifest: {output / 'manifest.json'}")
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=ROOT / "configs/validation_v1.yaml")
    args = parser.parse_args()
    run(args.config)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
