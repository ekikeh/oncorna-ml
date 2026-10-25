"""Pure generated-data checks of the reusable prediction input contract."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from oncorna.prediction import (
    FROZEN_INPUT_GENE_COUNT,
    PredictionInputContract,
    PredictionInputError,
)


@pytest.fixture
def frame() -> pd.DataFrame:
    return pd.DataFrame(
        [[0.0, 2.5, 3.0], [1.2, 0.5, 2.0]],
        index=["SYN-BRCA-0001", "SYN-BRCA-0002"],
        columns=["SYN-GENE-01", "SYN-GENE-02", "SYN-GENE-03"],
    )


def test_validation_returns_same_frame_without_transformation(frame: pd.DataFrame) -> None:
    before = frame.copy(deep=True)
    contract = PredictionInputContract(tuple(frame.columns))
    assert contract.validate(frame) is frame
    pd.testing.assert_frame_equal(frame, before)
    schema = list(frame.columns)
    from_list = PredictionInputContract(schema)
    schema.reverse()
    assert from_list.gene_ids == tuple(frame.columns)


@pytest.mark.parametrize("case", ["missing", "extra", "duplicate", "reordered", "transposed"])
def test_gene_schema_mismatch_fails(frame: pd.DataFrame, case: str) -> None:
    contract = PredictionInputContract(tuple(frame.columns))
    if case == "missing":
        frame = frame.iloc[:, :-1]
    elif case == "extra":
        frame["SYN-GENE-EXTRA"] = 1.0
    elif case == "duplicate":
        frame.columns = ["SYN-GENE-01", "SYN-GENE-01", "SYN-GENE-03"]
    elif case == "reordered":
        frame = frame.iloc[:, ::-1]
    else:
        frame = frame.T
    with pytest.raises(PredictionInputError, match="gene|Gene"):
        contract.validate(frame)


@pytest.mark.parametrize("sample_id", ["", " ", "SYN BAD", "SYN\nBAD", "../SYN", 1, None])
def test_invalid_sample_identifiers_fail(frame: pd.DataFrame, sample_id: object) -> None:
    frame.index = [sample_id, "SYN-BRCA-0002"]
    with pytest.raises(PredictionInputError, match="sample identifier"):
        PredictionInputContract(tuple(frame.columns)).validate(frame)


def test_duplicate_sample_identifiers_fail(frame: pd.DataFrame) -> None:
    frame.index = ["SYN-BRCA-0001"] * 2
    with pytest.raises(PredictionInputError, match="Duplicate sample"):
        PredictionInputContract(tuple(frame.columns)).validate(frame)


def test_multiindex_samples_fail(frame: pd.DataFrame) -> None:
    frame.index = pd.MultiIndex.from_tuples([("SYN", "1"), ("SYN", "2")])
    with pytest.raises(PredictionInputError, match="simple index"):
        PredictionInputContract(tuple(frame.columns)).validate(frame)


def test_unhashable_sample_and_gene_names_fail_cleanly(frame: pd.DataFrame) -> None:
    contract = PredictionInputContract(tuple(frame.columns))
    bad_samples = frame.copy()
    sample_names = np.empty(2, dtype=object)
    sample_names[:] = [["SYN", "1"], ["SYN", "2"]]
    bad_samples.index = pd.Index(sample_names)
    with pytest.raises(PredictionInputError, match="sample identifier"):
        contract.validate(bad_samples)
    gene_names = np.empty(3, dtype=object)
    gene_names[:] = [["SYN", "1"], ["SYN", "2"], ["SYN", "3"]]
    frame.columns = pd.Index(gene_names)
    with pytest.raises(PredictionInputError, match="gene identifiers"):
        contract.validate(frame)


@pytest.mark.parametrize("dtype", [str, bool, complex, object])
def test_nonnumeric_or_nonreal_dtypes_fail(frame: pd.DataFrame, dtype: type) -> None:
    invalid = frame.astype(dtype)
    with pytest.raises(PredictionInputError, match="real numeric"):
        PredictionInputContract(tuple(frame.columns)).validate(invalid)


@pytest.mark.parametrize("value", [np.nan, np.inf, -np.inf, -0.1])
def test_nonfinite_or_negative_values_fail(frame: pd.DataFrame, value: float) -> None:
    frame.iloc[0, 0] = value
    with pytest.raises(PredictionInputError, match="NaN|nonnegative"):
        PredictionInputContract(tuple(frame.columns)).validate(frame)


def test_nullable_missing_values_fail(frame: pd.DataFrame) -> None:
    frame = frame.astype("Float64")
    frame.iloc[0, 0] = pd.NA
    with pytest.raises(PredictionInputError, match="NaN"):
        PredictionInputContract(tuple(frame.columns)).validate(frame)


@pytest.mark.parametrize("shape", [(), (3,), (2, 3), (2, 3, 1)])
def test_arrays_without_explicit_axis_identifiers_fail(shape: tuple[int, ...]) -> None:
    with pytest.raises(PredictionInputError, match="two-dimensional.*DataFrame"):
        PredictionInputContract(("SYN-GENE-01",)).validate(np.zeros(shape))


def test_empty_sample_population_fails(frame: pd.DataFrame) -> None:
    with pytest.raises(PredictionInputError, match="at least one sample"):
        PredictionInputContract(tuple(frame.columns)).validate(frame.iloc[:0])


@pytest.mark.parametrize(
    "schema", [(), "SYN-GENE-01", {"SYN-GENE-01"}, ("A", "A"), ("",), (1,), ("A B",), ("A\x00",)]
)
def test_invalid_explicit_schema_fails(schema: object) -> None:
    with pytest.raises(PredictionInputError):
        PredictionInputContract(schema)


@pytest.mark.parametrize("case", ["missing", "extra", "duplicate", "reordered", "nan"])
def test_generated_frozen_size_contract_without_real_model(case: str) -> None:
    genes = tuple(f"SYN-SCHEMA-{index:05d}" for index in range(FROZEN_INPUT_GENE_COUNT))
    contract = PredictionInputContract.for_frozen_model(genes)
    frame = pd.DataFrame(np.ones((1, len(genes))), columns=genes, index=["SYN-BRCA-0001"])
    assert contract.validate(frame) is frame
    if case == "missing":
        frame = frame.iloc[:, :-1]
    elif case == "extra":
        frame["SYN-EXTRA"] = 1.0
    elif case == "duplicate":
        frame.columns = [genes[0], *genes[:-1]]
    elif case == "reordered":
        frame = frame.iloc[:, ::-1]
    else:
        frame.iloc[0, 0] = np.nan
    with pytest.raises(PredictionInputError):
        contract.validate(frame)


@pytest.mark.parametrize("count", [12, 17_613, 17_623, 20_529, 20_531])
def test_frozen_size_contract_rejects_wrong_universe(count: int) -> None:
    genes = tuple(f"SYN-SCHEMA-{index:05d}" for index in range(count))
    with pytest.raises(PredictionInputError, match="20,530"):
        PredictionInputContract.for_frozen_model(genes)
