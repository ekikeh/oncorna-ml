# Standalone synthetic demonstration

This teaching and software-validation fixture runs offline, creates all values in memory,
and prints aggregate JSON with the scope `synthetic_demo_only`. It does not read TCGA
files, the real training preview, or fitted artifacts. It does not import the protected
evaluation runner or authorization infrastructure. The demo writes its report to standard
output and does not save datasets or models.

## Run

From an installed checkout with Python 3.11 and the project's dependencies:

```bash
python -m oncorna.synthetic_demo
python -m oncorna.synthetic_demo --seed 20261026
```

The output includes the seed, ordered synthetic genes and classes, dataset digests,
dependency versions, aggregate synthetic accuracy and macro F1, and a 5-by-5 confusion
matrix. Rows are true synthetic classes and columns are predicted synthetic classes,
both in the stated class order. No sample-level predictions or risk estimates are printed.
Invalid seeds and validation failures exit with status 2 and a clear error message.

There are deliberately no input-path, output-path, model, or authorization options.
The only model is a fresh in-memory demonstration fit. A failed fit raises an error;
there is no model discovery, deserialization, download, or fallback model.

## Generated data and learning procedure

- Default seed: **20261026**, using NumPy `default_rng`.
- Samples: **100 training** (20 per class) and **40 evaluation** (8 per class).
- Invented sample IDs: `SYN-BRCA-0001` through `SYN-BRCA-0140`, with no overlap
  between the training and evaluation populations. Row order is shuffled within each.
- Ordered genes: `SYN-GENE-01` through `SYN-GENE-12`.
- Labels: `SYN-CLASS-A` through `SYN-CLASS-E`, not PAM50 calls.
- Values: finite, nonnegative, normalized-log-expression-like numbers. Baselines are
  drawn uniformly from 0.25 to 1.5; each class gets an added artificial signal in one
  of the first five columns and one of the next five columns. Two columns are noise.

The generator directly constructs expression-like values; these are not raw counts,
observed TCGA measurements, or a realistic simulation of breast-cancer biology. It uses
no cohort-derived parameters. A `StandardScaler` and L2 multinomial logistic regression
(`C=1`, `lbfgs`) fit only the 100 synthetic training samples. Evaluation only transforms
and predicts with that fitted pipeline. There is no hyperparameter search or calibration.

Repeated generation and execution with the same seed and dependency environment produce
the same values and predictions. JSON excludes timestamps and run durations and records
data digests and Python, NumPy, pandas, and scikit-learn versions. The separate
[packaging guide](packaging.md) documents exact platform locks and container/CI checks. Bitwise equality across
different library versions or numerical backends is not promised.

## Reusable input contract

`oncorna.prediction.PredictionInputContract` requires an explicit ordered gene list.
`validate(frame)` accepts only a nonempty two-dimensional pandas DataFrame with samples
as rows, gene identifiers as columns, and sample identifiers as the index. It returns
the original frame unchanged. It rejects missing, extra, duplicate, or reordered genes;
invalid or duplicate samples; NaN, infinity, negative values; unlabeled arrays; empty
populations; and hierarchical sample indices.

Supported column dtypes are native-endian NumPy `float32`, `float64`, and signed or
unsigned `int8`/`int16`/`int32`/`int64` and `uint8`/`uint16`/`uint32`/`uint64`.
Floating values must be finite and nonnegative. Integer values must lie in `0..2**53`
(inclusive), the contiguous integer range exactly representable as float64.
Finite float64 values up to its maximum and positive subnormal values are supported.
IEEE **negative zero is accepted as numerically zero** and its original sign bit is
preserved. Native values are checked column by column before any float64 conversion;
the conversion is checked for nonfinite or negative results and is only a temporary
validation copy. The returned frame retains its original values and dtypes.

Extended-precision floats, float16, non-native byte order, pandas nullable/sparse
extension dtypes, and complex, object, string, boolean, or structured dtypes are outside
this deliberately narrow contract. A platform alias with the same dtype as float64 is
treated as float64; tests do not claim extended-precision coverage where `longdouble`
has no wider exponent range. No dtype or string conversion is performed to admit
unsupported inputs. Numeric-looking strings are rejected rather than coerced.

Sample IDs must start with
an ASCII letter or digit and otherwise contain only ASCII letters, digits, underscores,
dots, colons, or hyphens. Schema genes must be unique nonempty printable strings without
whitespace. Validation never rearranges genes, normalizes, log-transforms, selects
features, fits a model, or verifies biological provenance.

```python
from oncorna.prediction import PredictionInputContract
from oncorna.synthetic_demo import DEMO_GENES, generate_synthetic_dataset

data = generate_synthetic_dataset()
checked = PredictionInputContract(DEMO_GENES).validate(data.evaluation_X)
```

`PredictionInputContract.for_frozen_model(ordered_genes)` additionally requires exactly
**20,530** genes, preserving the documented real-model input size. This API does not
supply or authenticate the real gene list. The caller must supply the appropriate
authoritative schema; having the right length does not prove identity or normalization.
Tests exercise this size with invented `SYN-SCHEMA-*` identifiers and generated values,
without constructing or loading the real model. The 12-gene demo schema and the
100-gene real training preview are not valid substitutes for that real input schema.

## Scientific and execution boundary

Synthetic scores demonstrate software behavior on deliberately easy planted patterns.
They are not research evidence, clinically calibrated risk, or an estimate of performance
on real expression. Do not compare them numerically with the completed TCGA-BRCA test
as equivalent experiments. The real study's 675-patient fit and closed 169-patient test
are described in the [model card](model_card.md) and [aggregate report](../reports/phase_3d_aggregate_results.md).

The actual fitted pipeline, full expression arrays, patient-level predictions, and
protected evaluation records are unavailable in the public package. The real evaluation
is never rerun by this demonstration. Windows guards for filesystem-backed synthetic
evaluation tests remain in force; this independent in-memory demo does not use those
protected APIs. See [reproducibility](reproducibility.md) and [release scope](release_scope.md).
