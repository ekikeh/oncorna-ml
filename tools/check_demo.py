"""Run two fresh demo processes and require complete same-environment JSON equality."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path


def check() -> dict:
    command = [sys.executable, "-B", "-m", "oncorna.synthetic_demo"]
    first = subprocess.check_output(command)
    second = subprocess.check_output(command)
    if first != second:
        raise ValueError("Complete synthetic demo output is not deterministic")
    result = json.loads(first)
    baseline = json.loads(
        (Path(__file__).resolve().parents[1] / "packaging/synthetic-baseline.json").read_text()
    )
    scientific_fields = {key: value for key, value in result.items() if key != "environment"}
    baseline_fields = {key: value for key, value in baseline.items() if key != "environment"}
    if scientific_fields != baseline_fields:
        raise ValueError(
            "Synthetic output changed from approved Phase 4C baseline; stop for review"
        )
    for dependency in ("numpy", "pandas", "scikit_learn"):
        if result["environment"][dependency] != baseline["environment"][dependency]:
            raise ValueError("Unreviewed numerical dependency version")
    if (
        result["scope"] != "synthetic_demo_only"
        or result["seed"] != 20261026
        or result["training_samples"] != 100
        or result["evaluation_samples"] != 40
        or len(result["gene_schema"]) != 12
        or len(result["class_order"]) != 5
        or sum(map(sum, result["confusion_matrix"]["counts"])) != 40
    ):
        raise ValueError("Unexpected synthetic demonstration contract")
    return result


if __name__ == "__main__":
    print(json.dumps(check(), sort_keys=True, indent=2, allow_nan=False))
