"""Validation only: samples-by-genes normalized, log-transformed expression.

No model loading, prediction, normalization, fitting, or gene rearrangement occurs here.
"""

from __future__ import annotations

import re
from dataclasses import dataclass

import numpy as np
import pandas as pd
from pandas.api.types import is_bool_dtype, is_complex_dtype, is_numeric_dtype

FROZEN_INPUT_GENE_COUNT = 20_530
INPUT_SEMANTICS = "previously_normalized_log_transformed_expression_like"
_SAMPLE_ID = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.:-]*", flags=re.ASCII)


class PredictionInputError(ValueError):
    """The supplied matrix or ordered schema violates the prediction contract."""


@dataclass(frozen=True)
class PredictionInputContract:
    """An explicit ordered gene universe; validation returns the original frame.

    The caller supplies gene identity. A size check alone cannot authenticate a real
    scientific schema or establish that upstream normalization was appropriate.
    """

    gene_ids: tuple[str, ...]

    def __post_init__(self) -> None:
        if not isinstance(self.gene_ids, (tuple, list)):
            raise PredictionInputError("Gene schema must be an explicit ordered tuple or list")
        genes = tuple(self.gene_ids)
        if not genes or any(
            not isinstance(gene, str)
            or not gene
            or any(char.isspace() for char in gene)
            or not gene.isprintable()
            for gene in genes
        ):
            raise PredictionInputError("Gene identifiers must be nonempty printable strings")
        if len(set(genes)) != len(genes):
            raise PredictionInputError("Gene schema contains duplicate gene identifiers")
        object.__setattr__(self, "gene_ids", genes)

    @classmethod
    def for_frozen_model(cls, gene_ids: tuple[str, ...]) -> PredictionInputContract:
        """Require 20,530 caller-supplied genes without opening any real artifact."""
        contract = cls(gene_ids)
        if len(contract.gene_ids) != FROZEN_INPUT_GENE_COUNT:
            raise PredictionInputError("Frozen model input requires exactly 20,530 ordered genes")
        return contract

    def validate(self, expression: pd.DataFrame) -> pd.DataFrame:
        """Reject malformed inputs; never repair, reorder, transform, or fit them."""
        if not isinstance(expression, pd.DataFrame) or expression.ndim != 2:
            raise PredictionInputError("Input must be a two-dimensional samples-by-genes DataFrame")
        if expression.empty or isinstance(expression.index, pd.MultiIndex):
            raise PredictionInputError("Input requires at least one sample with a simple index")
        if any(
            not isinstance(sample, str) or _SAMPLE_ID.fullmatch(sample) is None
            for sample in expression.index
        ):
            raise PredictionInputError(
                "Invalid sample identifier: use letters, digits, underscores, "
                "dots, colons or hyphens"
            )
        if expression.index.has_duplicates:
            raise PredictionInputError("Duplicate sample identifiers are not allowed")
        if any(not isinstance(gene, str) for gene in expression.columns):
            raise PredictionInputError("Invalid gene identifiers: column names must be strings")
        if expression.columns.has_duplicates:
            raise PredictionInputError("Duplicate gene identifiers are not allowed")
        if tuple(expression.columns) != self.gene_ids:
            raise PredictionInputError(
                "Gene columns must exactly match the supplied schema: "
                "missing, extra or reordered genes"
            )
        if any(
            not is_numeric_dtype(dtype) or is_bool_dtype(dtype) or is_complex_dtype(dtype)
            for dtype in expression.dtypes
        ):
            raise PredictionInputError("Expression values must have real numeric dtypes")
        values = expression.to_numpy(dtype=np.float64, na_value=np.nan)
        if not np.isfinite(values).all():
            raise PredictionInputError("NaN and infinite expression values are not allowed")
        if (values < 0).any():
            raise PredictionInputError("Expected nonnegative normalized log-expression-like values")
        return expression
