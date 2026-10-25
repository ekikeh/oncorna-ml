"""Standalone, in-memory teaching fixture; no real artifacts or evaluator imports."""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
from dataclasses import dataclass
from numbers import Integral
from typing import Any

import numpy as np
import pandas as pd
import sklearn
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, confusion_matrix, f1_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from oncorna.prediction import INPUT_SEMANTICS, PredictionInputContract

SCOPE = "synthetic_demo_only"
DEMO_SEED = 20261026
DEMO_GENES = tuple(f"SYN-GENE-{index:02d}" for index in range(1, 13))
DEMO_CLASSES = tuple(f"SYN-CLASS-{letter}" for letter in "ABCDE")
TRAIN_PER_CLASS = 20
EVALUATION_PER_CLASS = 8


@dataclass(frozen=True)
class SyntheticDataset:
    """Invented samples with disjoint training and evaluation identifiers."""

    train_X: pd.DataFrame
    train_y: pd.Series
    evaluation_X: pd.DataFrame
    evaluation_y: pd.Series
    seed: int


def generate_synthetic_dataset(seed: int = DEMO_SEED) -> SyntheticDataset:
    """Plant two artificial signals per class in nonnegative expression-like values.

    These are directly simulated log-scale-like numbers, not observed expression,
    a count model, or a realistic simulation of breast-cancer biology.
    """
    if isinstance(seed, bool) or not isinstance(seed, Integral) or not 0 <= seed < 2**32:
        raise ValueError("Seed must be an integer between 0 and 2**32 - 1")
    seed = int(seed)
    rng = np.random.default_rng(seed)

    def partition(per_class: int, first_id: int) -> tuple[pd.DataFrame, pd.Series]:
        classes = np.repeat(np.arange(len(DEMO_CLASSES)), per_class)
        values = rng.uniform(0.25, 1.5, (len(classes), len(DEMO_GENES)))
        for row, class_index in enumerate(classes):
            values[row, class_index] += rng.uniform(2.0, 3.0)
            values[row, class_index + len(DEMO_CLASSES)] += rng.uniform(1.0, 2.0)
        sample_ids = pd.Index(
            [f"SYN-BRCA-{index:04d}" for index in range(first_id, first_id + len(classes))],
            name="sample_id",
        )
        expression = pd.DataFrame(values, index=sample_ids, columns=DEMO_GENES)
        labels = pd.Series(
            [DEMO_CLASSES[index] for index in classes], index=sample_ids, name="synthetic_class"
        )
        order = rng.permutation(len(classes))
        return expression.iloc[order].copy(), labels.iloc[order].copy()

    train_X, train_y = partition(TRAIN_PER_CLASS, 1)
    evaluation_X, evaluation_y = partition(EVALUATION_PER_CLASS, len(train_X) + 1)
    return SyntheticDataset(train_X, train_y, evaluation_X, evaluation_y, seed)


def _data_digest(expression: pd.DataFrame, labels: pd.Series) -> str:
    header = json.dumps(
        {
            "samples": list(expression.index),
            "genes": list(expression.columns),
            "labels": list(labels),
        },
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    values = expression.to_numpy(dtype="<f8").tobytes(order="C")
    return hashlib.sha256(header + b"\n" + values).hexdigest()


def run_demo(seed: int = DEMO_SEED) -> dict[str, Any]:
    """Fit a fresh synthetic model and return aggregate results only.

    There is intentionally no model, artifact, input-path or output-path argument.
    No supplied or discovered model can substitute for this new demonstration fit.
    """
    data = generate_synthetic_dataset(seed)
    contract = PredictionInputContract(DEMO_GENES)
    train = contract.validate(data.train_X)
    evaluation = contract.validate(data.evaluation_X)
    pipeline = Pipeline(
        [
            ("scaler", StandardScaler()),
            (
                "classifier",
                LogisticRegression(C=1.0, solver="lbfgs", max_iter=500, random_state=data.seed),
            ),
        ]
    )
    pipeline.fit(train, data.train_y)
    # Prediction only transforms with the existing scaler; no fitting is done here.
    predicted = pipeline.predict(evaluation)
    return {
        "scope": SCOPE,
        "seed": data.seed,
        "input_semantics": INPUT_SEMANTICS,
        "gene_schema": list(DEMO_GENES),
        "class_order": list(DEMO_CLASSES),
        "training_samples": len(train),
        "evaluation_samples": len(evaluation),
        "training_data_sha256": _data_digest(train, data.train_y),
        "evaluation_data_sha256": _data_digest(evaluation, data.evaluation_y),
        "model": {
            "source": "fresh_synthetic_training_only",
            "steps": ["StandardScaler", "L2 multinomial LogisticRegression(C=1)"],
        },
        "metrics": {
            "macro_f1": float(
                f1_score(data.evaluation_y, predicted, labels=list(DEMO_CLASSES), average="macro")
            ),
            "accuracy": float(accuracy_score(data.evaluation_y, predicted)),
        },
        "confusion_matrix": {
            "rows": "true_synthetic_class",
            "columns": "predicted_synthetic_class",
            "class_order": list(DEMO_CLASSES),
            "counts": confusion_matrix(
                data.evaluation_y, predicted, labels=list(DEMO_CLASSES)
            ).tolist(),
        },
        "environment": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scikit_learn": sklearn.__version__,
        },
        "limitations": (
            "Teaching fixture only; synthetic metrics are not scientific or clinical results."
        ),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=f"Run the offline {SCOPE} teaching fixture")
    parser.add_argument("--seed", type=int, default=DEMO_SEED, help="Synthetic generator seed")
    args = parser.parse_args(argv)
    try:
        result = run_demo(args.seed)
    except ValueError as exc:
        parser.error(f"Synthetic demo input validation failed: {exc}")
    print(json.dumps(result, sort_keys=True, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
