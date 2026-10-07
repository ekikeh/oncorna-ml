# TCGA-BRCA classification cohort definition

**Generated:** 2026-10-07<br>
**Status:** Phase 1C cohort definition only. No expression QC, gene filtering, splitting, differential expression, GSEA, or machine learning was performed.

## Cohort rules

The primary outcome is the clinical field `PAM50Call_RNAseq`. The source value is copied unchanged to `pam50_original_label`; the separate normalized field maps `LumA` → `Luminal A`, `LumB` → `Luminal B`, `Basal` → `Basal-like`, `Her2` → `HER2-enriched`, and `Normal` → `Normal-like`. All five classes are retained.

Sample type is read from the fourth TCGA barcode field (the first two digits are the sample-type code) and checked against the clinical `sample_type_id` and `sample_type` fields. GDC defines code 01 as Primary Solid Tumor, 06 as Metastatic, and 11 as Solid Tissue Normal ([sample-type code table](https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/sample-type-codes); [TCGA barcode guide](https://docs.gdc.cancer.gov/Encyclopedia/pages/TCGA_Barcode/)). The downloaded clinical text labels for those codes are `Primary Tumor`, `Metastatic`, and `Solid Tissue Normal`.

The matrix also contains `bcr_sample_barcode`, `tumor_tissue_site`, `tissue_source_site`, `tissue_prospective_collection_indicator`, and `tissue_retrospective_collection_indicator`. These describe barcode, anatomic site, source site, or collection context; they do not replace the explicit sample-type code/text for distinguishing primary tumor, metastasis, and solid-tissue normal samples.

For each participant derived from the first three sample-barcode fields, the exact selection hierarchy is: **01 Primary Tumor → 06 Metastatic → 11 Solid Tissue Normal**. Within a type-priority tie, the lexicographically smallest full sample ID wins. No expression values, PAM50 labels/confidence, survival, or downstream results are used to choose a sample.

One sample is selected per patient before applying PAM50 availability. If the selected sample has missing PAM50, that patient is excluded; the builder does not fall back to a secondary sample with a nonmissing label. This preserves the primary-tumor preference. Unrecognized or conflicting sample-type evidence and conflicting `_PATIENT` values are flagged for human review rather than ranked automatically.

## Starting data and identifier checks

- Expression sample columns at start: **1,218**.
- Clinical sample rows: **1,247**; unique nonblank clinical `sampleID`s: **1,247**.
- Exact expression/clinical `sampleID` matches: **1,218 / 1,218**.
- Expression samples without an exact clinical match: **0**.
- Clinical-only sample IDs (not in expression): **29**; these are outside the expression-sample denominator.
- Duplicate clinical sample IDs: **0**.
- Unique expression participants from TCGA barcode parsing: **1,097**.
- Patients with multiple expression samples: **118**.

### Expression sample types

| Barcode code | Clinical `sample_type` | Matched expression samples | Missing PAM50 among these samples |
|---|---|---:|---:|
| `01` | `Primary Tumor` | 1,097 | 253 |
| `06` | `Metastatic` | 7 | 2 |
| `11` | `Solid Tissue Normal` | 114 | 7 |

For all matched expression IDs, the barcode code, clinical `sample_type_id`, and clinical `sample_type` text agree. Every derived participant ID also agrees with clinical `_PATIENT` in the matched rows. The older `patient_id` and `bcr_patient_barcode` fields are not used as the grouping key because they contain missing and truncated values (details below).

### TCGA identifier cross-checks on matched expression rows

| Clinical field/check | Matches | Mismatches | Missing |
|---|---:|---:|---:|
| Derived patient ID vs `_PATIENT` | 1,218 | 0 | 0 |
| Derived patient suffix vs `patient_id` | 1,213 | 2 | 3 |
| Derived patient ID vs `bcr_patient_barcode` | 1,213 | 2 | 3 |
| Derived patient and sample type vs `bcr_sample_barcode` | 1,024 | 2 | 192 |

Two expression rows (`TCGA-E9-A1NA-01` and `TCGA-E9-A1NA-11`) for participant `TCGA-E9-A1NA` have `patient_id`=`A1` and `bcr_patient_barcode`=`TCGA-E9-A1`, while `sampleID` and `_PATIENT` identify `TCGA-E9-A1NA`. The same two `bcr_sample_barcode` values use the truncated participant portion. `_PATIENT`, barcode sample-type code, `sample_type_id`, and `sample_type` agree for both rows. These rows are grouped using the sample ID's first three fields, cross-checked with `_PATIENT`; the truncated legacy fields are not used for grouping. Human review of this source inconsistency is still warranted.

- Sample IDs with conflicting legacy barcode fields: `TCGA-E9-A1NA-01`, `TCGA-E9-A1NA-11`.

Sample-type evidence status across matched expression rows: `consistent` 1,218.

## One-patient-one-sample resolution

- The expression sample set represents **1,097** participants; **118** have multiple samples.
- Samples per patient before selection: 1 sample(s): 979 patient(s); 2 sample(s): 115 patient(s); 3 sample(s): 3 patient(s).
- Participants with exactly one observed type-01 primary tumor sample: **1,097**.
- Type compositions among multi-sample participants:
  - `01+06`: 4 participant(s)
  - `01+06+11`: 3 participant(s)
  - `01+11`: 111 participant(s)

The 121 non-primary samples are not selected because every represented participant has a type-01 primary tumor sample. Of these 121 secondary samples, 112 have a nonmissing PAM50 label and 9 are missing it; neither label availability nor label identity changes the primary-tumor selection.

## Exclusions and final cohort

Among all **1,218** expression samples with matching clinical metadata, **262** have missing `PAM50Call_RNAseq` across all sample types. After the sample-type hierarchy selects one sample per patient, **253** selected primary-tumor samples are excluded for missing PAM50. The other **9** missing labels are on secondary samples already superseded by the primary-tumor rule and are not double-counted.

Expression-sample disposition (counts sum to the starting expression sample count):

| Reason | Count |
|---|---:|
| Retained: selected primary tumor with a mapped PAM50 label (`retained_primary_tumor`) | 844 |
| Retained: selected metastatic sample with a mapped PAM50 label (`retained_metastatic`) | 0 |
| Retained: selected solid-tissue-normal sample with a mapped PAM50 label (`retained_solid_tissue_normal`) | 0 |
| Excluded: selected preferred sample has missing PAM50 (`excluded_missing_pam50`) | 253 |
| Review: selected sample has a nonmissing, unmapped PAM50 label (`unmapped_pam50_label_requires_review`) | 0 |
| Not selected: a higher-priority sample type exists for this patient (`not_selected_lower_priority_sample_type`) | 121 |
| Not selected: lexicographically later ID among same-priority samples (`not_selected_lexicographic_tie_break`) | 0 |
| Excluded: no exact clinical sampleID match (`no_matching_clinical_metadata`) | 0 |
| Review: clinical sampleID occurs more than once (`duplicate_clinical_sample_id_requires_review`) | 0 |
| Review: patient ID could not be derived from TCGA barcode (`unparseable_tcga_patient_id`) | 0 |
| Review: clinical _PATIENT conflicts with sample barcode (`patient_identifier_conflict_requires_review`) | 0 |
| Review: sample-type evidence is conflicting or unmapped (`sample_type_conflict_requires_review`) | 0 |

- Clinical-only sample IDs outside the expression denominator: **29**.

- Selected one-sample-per-patient records before PAM50 eligibility: **1,097**.
- Final classification cohort: **844 unique patients**.

### Final normalized PAM50 distribution

| Normalized PAM50 label | Samples/patients |
|---|---:|
| `Luminal A` | 421 |
| `Luminal B` | 192 |
| `Basal-like` | 141 |
| `HER2-enriched` | 67 |
| `Normal-like` | 23 |

All five classes are retained, including `Normal-like`; no minimum class-size rule has been applied. `classification_cohort.tsv` contains metadata only. `cohort_sample_audit.tsv` contains one disposition row per expression sample so every selection and exclusion is traceable. Both files are written under the Git-ignored `data/processed/` directory; the full expression matrix is not copied.

## Reproducibility and scope boundary

Regenerate the cohort tables and this document from the repository root with `python scripts/build_cohort.py`. The builder reads the expression header only, matches clinical `sampleID`s exactly after outer-whitespace trimming, and writes metadata only. It does not perform expression QC, gene filtering, train/test splitting, differential expression, GSEA, or machine learning.

## Human review items

1. Review the truncated `TCGA-E9-A1NA` legacy clinical barcode fields described above; the grouping key is supported by `_PATIENT` and expression `sampleID`.
2. Confirm that selecting the primary tumor before PAM50 filtering, with no fallback when that selected label is missing, remains the desired rule for later work.
3. Confirm any future sample-type codes outside 01/06/11 before assigning them a selection rank; this build does not guess an order for unrecognized codes.
