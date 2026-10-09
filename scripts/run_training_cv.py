"""Run the approved Phase 2B training-only CV; never load held-out expression rows."""

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

import numpy as np
import pandas as pd
import sklearn
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from oncorna.modeling import (  # noqa: E402
    evaluate_candidate,
    load_training_data,
    make_folds,
    select_c,
)
from oncorna.preprocessing import sha256_file  # noqa: E402


def write_json(path: Path, value: Any) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(value, handle, indent=2, sort_keys=True, ensure_ascii=False)
        handle.write("\n")


def git_value(*args: str) -> str:
    return subprocess.check_output(["git", *args], cwd=ROOT, text=True).strip()


def run(config_path: Path) -> dict[str, Any]:
    config_path = config_path.resolve()
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    if (
        config["version"] != "modeling_v1"
        or config["cv_folds"] != 4
        or config["cv_seed"] != 20261009
        or config["c_values"] != [0.001, 0.01, 0.1, 1.0]
        or config["solver"] != "lbfgs"
        or config["penalty"] != "l2"
        or config["class_weight"] is not None
        or config["n_jobs"] != 1
        or config["primary_metric"] != "macro_f1"
        or config["zero_division"] != 0
        or config["labels"]
        != ["Luminal A", "Luminal B", "Basal-like", "HER2-enriched", "Normal-like"]
        or config["expression_gt"] != 1.0
        or config["minimum_sample_fraction"] != 0.20
        or config["training_size"] != 506
        or config["gene_count"] != 20530
    ):
        raise ValueError("Configuration differs from approved Phase 2B protocol")
    output = (ROOT / config["output_dir"]).resolve()
    processed = (ROOT / "data/processed").resolve()
    if not output.is_relative_to(processed) or output == processed or output.exists():
        raise ValueError("Output must be a new directory under local-only data/processed")
    data = load_training_data(ROOT, config)
    folds = make_folds(data, config["cv_folds"], config["cv_seed"])
    output.mkdir(parents=True)
    assignments = {
        "seed": config["cv_seed"],
        "n_folds": config["cv_folds"],
        "folds": folds,
    }
    write_json(output / "fold_assignments.json", assignments)
    results = []
    predictions = []
    started_at = datetime.now(timezone.utc).isoformat()
    for c_value in [None, *config["c_values"]]:
        name = "dummy" if c_value is None else f"logistic_C_{c_value:g}"
        print(f"Running {name} (4 folds)", flush=True)
        try:
            result, rows = evaluate_candidate(data, folds, config, c_value)
        except Exception as exc:
            write_json(
                output / "failure.json",
                {
                    "candidate": name,
                    "error_type": type(exc).__name__,
                    "message": str(exc),
                    "traceback": traceback.format_exc(),
                    "completed_candidates": [item["C"] for item in results],
                },
            )
            raise
        results.append(result)
        predictions.extend(rows)
        write_json(output / "candidate_results.json", results)
        print(f"  mean fold macro F1: {result['fold_macro_f1_mean']:.4f}", flush=True)
    selected_c = select_c(results)
    with (output / "oof_predictions.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            fieldnames=[
                "patient_id",
                "sample_id",
                "true_label",
                "predicted_label",
                "model",
                "C",
            ],
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(predictions)
    manifest = {
        "version": "modeling_v1",
        "scope": "training_only_cv",
        "started_at_utc": started_at,
        "finished_at_utc": datetime.now(timezone.utc).isoformat(),
        "git_commit_at_run": git_value("rev-parse", "HEAD"),
        "git_dirty_at_run": bool(git_value("status", "--porcelain")),
        "config_path": config_path.relative_to(ROOT).as_posix(),
        "config_sha256": sha256_file(config_path),
        "input_hashes": data.input_hashes,
        "training_count": len(data.y),
        "gene_count": data.X.shape[1],
        "class_order": config["labels"],
        "class_counts": {label: int(sum(data.y == label)) for label in config["labels"]},
        "cv_seed": config["cv_seed"],
        "cv_folds": config["cv_folds"],
        "fold_assignments_sha256": sha256_file(output / "fold_assignments.json"),
        "selected_C": selected_c,
        "selection_rule": "maximum mean fold macro F1; exact tie goes to smaller C",
        "selection_bias_note": "CV selected C; development metrics may be optimistic",
        "software": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scikit_learn": sklearn.__version__,
        },
        "output_hashes": {
            name: sha256_file(output / name)
            for name in ("fold_assignments.json", "candidate_results.json", "oof_predictions.tsv")
        },
    }
    write_json(output / "manifest.json", manifest)
    print(f"Selected C: {selected_c:g}", flush=True)
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=ROOT / "configs/modeling_v1.yaml")
    args = parser.parse_args()
    run(args.config)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
