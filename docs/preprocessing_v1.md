# Training-only gene preprocessing, version 1

This foundation corrects the Phase 1D/1E cohort-wide gene-prevalence screen for
future PAM50 classification. It does not train a model or change the approved
cohort, frozen `split_v1`, historical 17,623-gene matrix, QC outputs, or demo.

## Inputs and boundary

`scripts/prepare_preprocessing_v1.py` reads the original
`data/raw/TCGA-BRCA_HiSeqV2.tsv.gz` and selects the 844 approved sample IDs
from `data/processed/classification_cohort.tsv`. It retains all 20,530 source
gene rows in source order. The source values are already
`log2(normalized_count + 1)`; the command parses them as float64 without
another expression transform. It checks the raw source against the tracked
`data/MANIFEST.sha256` and the versioned config's expected checksum.

The existing `data/processed/split_v1.json` is read and SHA-256 checked before
use. The command validates each recorded patient/sample pair, partition count,
cohort order, class count, completeness, and disjointness. It also checks the
cohort and historical Phase 1E matrix digests recorded in the frozen manifest.
It never calls the split generator. Validation and test expression values are
copied into the unfiltered local view but are not used to learn gene retention.

## Fixed gene filter

`GenePrevalenceFilter.fit(X)` requires a numeric pandas DataFrame with samples
as rows and unique gene IDs as columns. It ignores `y`. A gene is retained if
its finite expression is **strictly greater than 1** in at least
`ceil(0.20 × n_fit)` fitting samples. Missing values do not pass the threshold;
the denominator remains all fitting samples. Infinite values are rejected.
With the 506 frozen training patients, the required count is 102. Gene order
remains source order, and the fitted mask is applied unchanged by `transform`
to validation and test. A mismatched input gene schema fails explicitly.

For later cross-validation, place this transformer inside a scikit-learn
`Pipeline` and fit a fresh pipeline on each fold's training patients. The
required count then depends on that fold's fitting size. Never fit the filter
once on all 506 patients before cross-validating those same patients. The
frozen validation and test partitions must not enter any fold's fit. This
checkpoint adds no imputer, scaler, feature selector, PCA, classifier, or
performance evaluation.

The numerical threshold is inherited from Phase 1 QC, where full-cohort
expression summaries had been inspected. Its historical selection context
should remain disclosed. This version keeps the approved rule fixed and does
not tune it from validation or test data.

## Local artifacts and command

From the repository root, with the Python 3.11 environment active:

```powershell
python scripts/prepare_preprocessing_v1.py
```

The command publishes a new `data/processed/preprocessing_v1/` directory:

- `all_gene_expression.npy`: float64 samples-by-genes values. Its expected
  payload is about 132 MiB; NumPy can reopen it with memory mapping on Windows.
- `all_gene_schema.json`: ordered sample and gene identifiers, dimensions,
  orientation, units, and dtype.
- `training_gene_list.tsv`: ordered retained gene IDs and training prevalence
  counts.
- `manifest.json`: input and output SHA-256 digests, source and cohort
  provenance, frozen split digest, partition identity digests, training sample
  digest, filter settings and counts, config compatibility, transform-only
  checks, and Python/library versions.

All four files are local-only under the existing `data/processed/*` ignore
rule. The command refuses to overwrite an existing versioned output directory;
review and version any later method changes separately. It does not require
the new retained-gene count to equal Phase 1E's 17,623.

## Historical provenance on Windows

The frozen split records a raw SHA-256 of the Phase 1 configuration written
with LF line endings. Windows Git checked out the same content with CRLF
line endings. The read-only compatibility check accepts an exact raw match or
a proven LF/CRLF-only variant: it hashes the current raw bytes, computes the
line-ending variant, and requires that variant's hash to equal the frozen
record. It then compares parsed split-relevant settings. The frozen raw digest
is retained and the compatibility mode is recorded in the new manifest.
Changes beyond line endings fail.

The split manifest records scikit-learn 1.6.1 for its historical generation;
the current preprocessing environment may use a later compatible release.
This command validates the existing artifact and assignments without rerunning
split generation. Exact historical reproduction would require the original
configuration bytes and recorded software environment in a separate working
copy.

## Limitations

The Xena HiSeqV2 resource is already normalized upstream; this workflow
cannot independently establish how that upstream normalization used samples.
The approved cohort includes one sample per patient and retains two previously
flagged low-expression outliers. No outlier is removed here. PAM50 labels are
used only to validate the frozen cohort/split accounting, never to fit the
gene-prevalence filter. The output remains research data and is not a clinical
diagnostic.
