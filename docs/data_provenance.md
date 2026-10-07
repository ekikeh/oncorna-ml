# Data provenance — Phase 1A acquisition and Phase 1B audit

## Selected dataset

- **Dataset:** TCGA Breast Cancer (TCGA-BRCA), legacy TCGA cohort hosted by the UCSC Xena TCGA hub.
- **Public source:** [UCSC Xena TCGA hub](https://tcga.xenahubs.net/) and [Xena Data Pages](https://xenabrowser.net/datapages/).
- **Access level:** public/open access; the selected files download without login or an API token. This workflow does not request controlled-access data and does not download FASTQ or BAM files.
- **Download date:** 2026-10-06 (initial retrieval run in this workspace). File-level SHA-256 checksums are in `data/MANIFEST.sha256`. Update this date if documenting a later retrieval/release.
- **Source versions:** the Xena data pages identify the expression matrix as version `2017-10-13` and the clinical matrix as version `2019-12-06`. Both objects returned HTTP `Last-Modified: 2021-04-08` during the 2026-10-06 source check.

## Files and data types

| Local file | Xena dataset / source URL | Data type and notes |
| --- | --- | --- |
| `data/raw/TCGA-BRCA_HiSeqV2.tsv.gz` | [`TCGA.BRCA.sampleMap/HiSeqV2`](https://tcga.xenahubs.net/download/TCGA.BRCA.sampleMap/HiSeqV2.gz) · [dataset metadata](https://xenabrowser.net/datapages/?dataset=TCGA.BRCA.sampleMap%2FHiSeqV2&host=https%3A%2F%2Ftcga.xenahubs.net) | Gene-level Illumina HiSeq RNA-seq matrix. Xena page metadata lists 20,531 identifiers × 1,218 samples; the Phase 1B scan found 20,530 data rows × 1,218 sample columns, a one-identifier discrepancy that remains unresolved. Unit is `log2(norm_count + 1)`; these are **not raw integer counts**. |
| `data/raw/TCGA-BRCA_BRCA_clinicalMatrix.tsv` | [`TCGA.BRCA.sampleMap/BRCA_clinicalMatrix`](https://tcga.xenahubs.net/download/TCGA.BRCA.sampleMap/BRCA_clinicalMatrix) · [dataset metadata](https://xenabrowser.net/datapages/?dataset=TCGA.BRCA.sampleMap%2FBRCA_clinicalMatrix&host=https%3A%2F%2Ftcga.xenahubs.net) | Tab-separated clinical/sample annotations: 1,247 samples × 194 fields in the Xena page metadata. Includes `sampleID` and the `PAM50Call_RNAseq` field where available; the source label is not inferred by this project. |

The expression matrix is a processed legacy TCGA RNASeqV2 dataset, not GDC-harmonized raw counts. Its values are `log2(normalized_count + 1)`: it may be used for visualization or documented ML preprocessing, but **must not be passed directly to DESeq2/PyDESeq2 as raw integer counts**. If a later phase requires count-based differential expression, select and document a raw-count source separately rather than reversing this log transform and treating the result as raw counts.

The Xena pages report different sample counts for the expression and clinical tables. Phase 1A deliberately downloads them without matching, filtering, or preprocessing. Later work must explicitly inspect barcode overlap, tumor/normal sample types, duplicate aliquots, and missing PAM50 labels before analysis.

No filtering, sample matching, duplicate handling, or expression transformation is performed by the downloader. The clinical matrix contains TCGA barcodes and clinical fields; keep the raw clinical matrix local and never commit it to Git.

## Phase 1B audit and identifier harmonization

The read-only audit is reproducible from the repository root:

```bash
python scripts/audit_data.py
```

It streams the gzipped expression matrix row by row and writes the preliminary findings to `docs/data_audit.md`. On the files downloaded for this run, the local scan found 20,530 expression data rows and 1,218 sample columns, plus 1,247 clinical sample rows and 194 fields. It found 20,530 unique nonblank gene IDs (20,501 gene-symbol-like strings and 29 Xena numeric fallbacks; no duplicate gene IDs). No Ensembl-style IDs were detected by the visible-pattern check. The Xena page metadata lists 20,531 expression identifiers, one more than the local data-row count; the cause is unresolved, so the report uses the scanned file count.

The exact sample-ID intersection is 1,218: 100.0% of the 1,218 unique expression IDs and 97.7% of the 1,247 unique clinical IDs. There are no expression-only IDs and 29 clinical-only IDs; the complete unmatched list is in `docs/data_audit.md`. Participant IDs from expression barcodes are derived from the first three hyphen-separated fields without collapsing samples. The scan found 1,097 expression participants, of whom 118 have multiple expression samples (979 participants with one sample, 115 with two, and 3 with three). The clinical `patient_id` field is an abbreviated final barcode segment in all 1,242 rows where both it and `bcr_patient_barcode` are present; it is not a full TCGA patient barcode.

`PAM50Call_RNAseq` has 956 nonmissing calls and 291 missing values. Its original labels and counts are preserved in the audit report, along with possible display-name correspondences and the distinct `PAM50_mRNA_nature2012` field. No labels are recoded, genes or samples filtered, unmatched IDs discarded, or duplicate aliquots selected. The audit does not make a clinical endpoint or sample-retention decision.

## Usage, license, and attribution

The Xena data pages do not display a file-level license for these cohort files. Cite UCSC Xena and the applicable TCGA source/publications, and check current TCGA publication/data-use guidance before reusing or redistributing source data. The repository's MIT license applies to project code only; it does not relicense these cohort files. This project commits checksums and provenance, not the full expression or clinical matrices.

Useful source pages:

- [UCSC Xena Data Pages](https://xenabrowser.net/datapages/)
- [TCGA / GDC Data Portal](https://portal.gdc.cancer.gov/)
- [Xena platform information](https://xenabrowser.net/)

## Why this source was selected

1. Expression and clinical/PAM50 annotations are available as two direct files from the same public TCGA Xena hub, simplifying beginner-friendly, reproducible acquisition.
2. The clinical matrix contains a PAM50-from-RNA-seq field (`PAM50Call_RNAseq`) that can be inspected and matched to expression sample barcodes in a later phase.
3. The files are public processed tables and require no credentials; the workflow avoids controlled-access data, FASTQ, and BAM downloads.
4. Direct file URLs and checksums make it straightforward to detect source-file changes.

The tradeoff is important: `HiSeqV2` is an older normalized-expression matrix, **not a raw-count matrix**. It was chosen for this acquisition-only phase because it is a stable, easy-to-download expression source paired with the requested subtype annotations. This decision does not change the scientific design, and no modeling or differential-expression analysis is included here.

## Download and checksum workflow

From the repository root:

```bash
python scripts/download_data.py
```

The script creates `data/raw/`, streams the two public files to disk, computes SHA-256 digests, and writes/updates `data/MANIFEST.sha256`. Existing files are verified and skipped by default. To replace them:

```bash
python scripts/download_data.py --force
```

Full files under `data/raw/` and `data/processed/` are excluded by `.gitignore`; `.gitkeep` placeholders and `data/MANIFEST.sha256` remain trackable. Do not force-add cohort files to Git.
