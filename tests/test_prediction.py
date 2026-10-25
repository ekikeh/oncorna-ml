"""Pure generated-data checks of the reusable prediction input contract."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from oncorna.prediction import (
    FROZEN_INPUT_GENE_COUNT,
    MAX_EXACT_INTEGER,
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
    with pytest.raises(PredictionInputError, match="Unsupported expression dtype"):
        PredictionInputContract(tuple(frame.columns)).validate(invalid)


@pytest.mark.parametrize("value", [np.nan, np.inf, -np.inf, -0.1])
def test_nonfinite_or_negative_values_fail(frame: pd.DataFrame, value: float) -> None:
    frame.iloc[0, 0] = value
    with pytest.raises(PredictionInputError, match="NaN|nonnegative"):
        PredictionInputContract(tuple(frame.columns)).validate(frame)


def test_nullable_missing_values_have_unsupported_dtype(frame: pd.DataFrame) -> None:
    frame = frame.astype("Float64")
    frame.iloc[0, 0] = pd.NA
    with pytest.raises(PredictionInputError, match="Unsupported expression dtype"):
        PredictionInputContract(tuple(frame.columns)).validate(frame)


@pytest.mark.parametrize(
    "dtype",
    [
        "float32",
        "float64",
        "int8",
        "int16",
        "int32",
        "int64",
        "uint8",
        "uint16",
        "uint32",
        "uint64",
    ],
)
def test_allowlisted_native_dtypes_preserve_values_and_input(dtype: str) -> None:
    native = np.dtype(dtype)
    if native.kind == "f":
        values = [0.0, np.nextafter(native.type(0), native.type(1)), 0.125, np.finfo(native).max]
    else:
        values = [0, 1, min(np.iinfo(native).max, MAX_EXACT_INTEGER)]
    expression = pd.DataFrame(
        np.asarray(values, dtype=native).reshape(-1, 1),
        columns=["SYN-GENE-01"],
        index=[f"SYN-{index:04d}" for index in range(len(values))],
    )
    before = expression.copy(deep=True)
    assert PredictionInputContract(("SYN-GENE-01",)).validate(expression) is expression
    pd.testing.assert_frame_equal(expression, before)
    np.testing.assert_array_equal(
        expression.to_numpy(dtype=np.float64).astype(native), expression.to_numpy()
    )


@pytest.mark.parametrize("dtype", ["float32", "float64"])
@pytest.mark.parametrize("case", ["tiny_negative", "large_negative", "nan", "inf", "negative_inf"])
def test_native_float_edge_cases_are_rejected_without_mutation(dtype: str, case: str) -> None:
    native = np.dtype(dtype).type
    values = {
        "tiny_negative": -np.nextafter(native(0), native(1)),
        "large_negative": -np.finfo(native).max,
        "nan": np.nan,
        "inf": np.inf,
        "negative_inf": -np.inf,
    }
    expression = pd.DataFrame(
        np.array([[values[case]]], dtype=dtype), columns=["SYN-GENE-01"], index=["SYN-0001"]
    )
    before = expression.copy(deep=True)
    with pytest.raises(PredictionInputError, match="NaN|nonnegative"):
        PredictionInputContract(("SYN-GENE-01",)).validate(expression)
    pd.testing.assert_frame_equal(expression, before)


@pytest.mark.parametrize("dtype", ["float32", "float64"])
def test_negative_zero_is_accepted_as_zero_and_preserved(dtype: str) -> None:
    expression = pd.DataFrame(
        np.array([[-0.0], [0.0]], dtype=dtype),
        columns=["SYN-GENE-01"],
        index=["SYN-0001", "SYN-0002"],
    )
    before = expression.copy(deep=True)
    assert PredictionInputContract(("SYN-GENE-01",)).validate(expression) is expression
    assert np.signbit(expression.iloc[0, 0]) and not np.signbit(expression.iloc[1, 0])
    pd.testing.assert_frame_equal(expression, before)


@pytest.mark.parametrize("dtype", ["int64", "uint64"])
@pytest.mark.parametrize("case", ["rounding", "maximum"])
def test_integer_values_outside_exact_conversion_range_fail(dtype: str, case: str) -> None:
    value = MAX_EXACT_INTEGER + 1 if case == "rounding" else np.iinfo(dtype).max
    expression = pd.DataFrame(
        np.array([[value]], dtype=dtype), columns=["SYN-GENE-01"], index=["SYN-0001"]
    )
    before = expression.copy(deep=True)
    with pytest.raises(PredictionInputError, match=r"between 0 and 2\*\*53"):
        PredictionInputContract(("SYN-GENE-01",)).validate(expression)
    pd.testing.assert_frame_equal(expression, before)


@pytest.mark.parametrize("dtype", ["int8", "int16", "int32", "int64"])
def test_negative_integer_values_fail(dtype: str) -> None:
    expression = pd.DataFrame(
        np.array([[np.iinfo(dtype).min]], dtype=dtype),
        columns=["SYN-GENE-01"],
        index=["SYN-0001"],
    )
    with pytest.raises(PredictionInputError, match="nonnegative"):
        PredictionInputContract(("SYN-GENE-01",)).validate(expression)


@pytest.mark.parametrize("case", ["negative", "nonfinite", "mixed_integer_rounding"])
def test_native_checks_run_before_any_matrix_conversion(monkeypatch, case: str) -> None:
    expression = pd.DataFrame(
        {
            "SYN-GENE-01": np.array([1], dtype="int64"),
            "SYN-GENE-02": np.array([MAX_EXACT_INTEGER + 1], dtype="uint64"),
        },
        index=["SYN-0001"],
    )
    if case != "mixed_integer_rounding":
        expression = expression.astype("float64")
        expression.iloc[0, 0] = -np.nextafter(0.0, 1.0) if case == "negative" else np.inf
    before = expression.copy(deep=True)

    def forbidden_cast(*_args, **_kwargs):
        pytest.fail("Matrix conversion reached before invalid native values were rejected")

    with monkeypatch.context() as guard:
        guard.setattr(pd.DataFrame, "to_numpy", forbidden_cast)
        with pytest.raises(PredictionInputError, match="nonnegative|NaN|between 0"):
            PredictionInputContract(tuple(expression.columns)).validate(expression)
    pd.testing.assert_frame_equal(expression, before)


@pytest.mark.parametrize(
    "dtype", ["float16", "Float32", "Float64", "Int64", "UInt64", pd.SparseDtype("float64", 0)]
)
def test_other_numeric_dtypes_are_explicitly_unsupported(dtype: object) -> None:
    expression = pd.DataFrame({"SYN-GENE-01": pd.Series([1], dtype=dtype)})
    expression.index = ["SYN-0001"]
    with pytest.raises(PredictionInputError, match="Unsupported expression dtype"):
        PredictionInputContract(("SYN-GENE-01",)).validate(expression)


def test_non_native_byte_order_is_unsupported() -> None:
    dtype = np.dtype("float64").newbyteorder("S")
    expression = pd.DataFrame(
        np.array([[1.0]], dtype=dtype), columns=["SYN-GENE-01"], index=["SYN-0001"]
    )
    with pytest.raises(PredictionInputError, match="Unsupported expression dtype"):
        PredictionInputContract(("SYN-GENE-01",)).validate(expression)


def test_structured_arrays_are_rejected() -> None:
    values = np.array([(1.0,)], dtype=[("SYN-GENE-01", "float64")])
    with pytest.raises(PredictionInputError, match="DataFrame"):
        PredictionInputContract(("SYN-GENE-01",)).validate(values)


@pytest.mark.skipif(
    np.finfo(np.longdouble).max <= np.finfo(np.float64).max,
    reason="Platform longdouble has no extended exponent range; no underflow/overflow coverage",
)
@pytest.mark.parametrize("value", ["-1e-4000", "1e4000"])
def test_extended_precision_is_rejected_before_lossy_conversion(monkeypatch, value: str) -> None:
    native = np.longdouble(value)
    assert np.isfinite(native)
    with np.errstate(over="ignore", under="ignore"):
        narrowed = np.float64(native)
    if native < 0:
        assert narrowed == 0 and np.signbit(narrowed)
    else:
        assert np.isinf(narrowed)
    expression = pd.DataFrame(
        np.array([[native]], dtype=np.longdouble), columns=["SYN-GENE-01"], index=["SYN-0001"]
    )
    before = expression.copy(deep=True)

    def forbidden_cast(*_args, **_kwargs):
        pytest.fail("Unsupported extended precision reached a float64 conversion")

    with monkeypatch.context() as guard:
        guard.setattr(pd.DataFrame, "to_numpy", forbidden_cast)
        with pytest.raises(PredictionInputError, match="Unsupported expression dtype"):
            PredictionInputContract(("SYN-GENE-01",)).validate(expression)
    pd.testing.assert_frame_equal(expression, before)


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
