# TCGA-BRCA data audit and identifier harmonization

**Generated:** 2026-10-07<br>
**Scope:** Phase 1B data audit only. No samples or genes were filtered, no labels were changed, and no modeling/DE/GSEA was performed.

## Inputs and dimensions

| Input | Rows/items | Columns/items | Notes |
|---|---:|---:|---|
| Expression (`TCGA-BRCA_HiSeqV2.tsv.gz`) | 20,530 gene rows | 1,218 sample columns | First field: `sample` |
| Clinical (`TCGA-BRCA_BRCA_clinicalMatrix.tsv`) | 1,247 sample rows | 194 fields | Tab-separated; sample ID field: `sampleID` |

### Expression matrix structure

- Gene rows: **20,530**; nonblank identifiers: **20,530**; unique gene IDs: **20,530**.
- Expression sample columns: **1,218**; unique nonblank sample IDs: **1,218**; blank sample headers: **0**.
- First gene identifiers: `ARHGEF10L`, `HIF3A`, `RNF17`, `RNF10`, `RNF11`.
- First sample IDs: `TCGA-AR-A5QQ-01`, `TCGA-D8-A1JA-01`, `TCGA-BH-A0BQ-01`, `TCGA-BH-A0BT-01`, `TCGA-A8-A06X-01`.
- Gene identifier pattern counts (classifying visible strings only; no mapping applied):
  - Gene-symbol-like identifier: 20,501
  - Xena fallback with numeric ID: 29
- Ensembl-style identifiers detected by this pattern check: **0**.
- Duplicate gene identifiers: **none detected**.
- Duplicate exact expression sample IDs: **none detected**.
- Rows with a width different from the header: **0**.
- Xena Data Pages metadata lists 20,531 identifiers, while the downloaded file contains 20,530 data rows after its header (difference: 1). This report uses the scanned file row count; the cause of the discrepancy is unresolved. [Metadata](https://xenabrowser.net/datapages/?dataset=TCGA.BRCA.sampleMap%2FHiSeqV2&host=https%3A%2F%2Ftcga.xenahubs.net).

### Identifier systems

- Expression columns use TCGA sample barcodes (for example, `TCGA-AR-A5QQ-01`).
- Expression row identifiers are predominantly gene-symbol-like; any Ensembl-like or Xena fallback IDs are counted above. This is a pattern audit, not annotation remapping.
- The clinical matrix exposes `sampleID`, `patient_id`, `bcr_patient_barcode`, and `_PATIENT`. The clinical `patient_id` field is not the full TCGA participant barcode: where both fields are present, it matches the final barcode segment. The full barcode fields retain the `TCGA-XX-XXXX` form.
- For expression barcodes, the participant-level ID uses the first three hyphen-separated fields joined unchanged: `TCGA-AR-A5QQ-01` → `TCGA-AR-A5QQ`. This keeps sample/aliquot suffixes out of the patient ID. No samples are collapsed by this audit.
- Comparison check for clinical `patient_id` suffix versus `bcr_patient_barcode`: 1,242 matches among 1,242 rows with both fields present.
- Example clinical mapping: `sampleID` `TCGA-3C-AAAU-01`, `patient_id` `AAAU` matches the final segment of `bcr_patient_barcode` `TCGA-3C-AAAU`.

## Expression-to-clinical sample matching

Sample IDs are compared exactly after trimming outer whitespace only. No barcode truncation, case conversion, or sample-type filtering is applied.

| Measure | Count |
|---|---:|
| Expression sample columns | 1,218 |
| Unique nonblank expression sample IDs | 1,218 |
| Clinical sample rows | 1,247 |
| Unique nonblank clinical sample IDs | 1,247 |
| Duplicate exact clinical sample IDs | 0 |
| Clinical rows with malformed width | 0 |
| Exact sample-ID intersection | 1,218 |
| Expression-only unique IDs | 0 |
| Clinical-only unique IDs | 29 |
| Expression IDs matched (intersection / unique expression IDs) | 100.0% |
| Clinical IDs matched (intersection / unique clinical IDs) | 97.7% |

Unmatched IDs are reported, not discarded:

- Expression-only: None.
- Clinical-only: `TCGA-A7-A4SA-11`, `TCGA-A7-A4SC-11`, `TCGA-A7-A4SD-11`, `TCGA-A7-A4SE-11`, `TCGA-AC-A2FG-11`, `TCGA-AC-A5EI-01`, `TCGA-AN-A0FE-01`, `TCGA-AN-A0FG-01`, `TCGA-BH-A0BF-11`, `TCGA-BH-A0BL-11`, `TCGA-BH-A0BO-11`, `TCGA-BH-A0DE-11`, `TCGA-BH-A0DI-11`, `TCGA-BH-A18F-11`, `TCGA-BH-A1ES-11`, `TCGA-BH-A1EY-11`, `TCGA-BH-A1F5-11`, `TCGA-C8-A9FZ-01`, `TCGA-E2-A14Y-11`, `TCGA-E2-A15L-11`, `TCGA-E2-A1B5-11`, `TCGA-E2-A1IF-11`, `TCGA-E2-A1II-11`, `TCGA-E2-A1IO-11`, `TCGA-E2-A1LI-11`, `TCGA-E9-A1N8-11`, `TCGA-E9-A1NE-11`, `TCGA-E9-A1NH-11`, `TCGA-XX-A899-11`.
- Clinical rows with missing sample IDs: **0**.

## Patient/sample relationship in expression

- Unique patients derived from expression sample barcodes: **1,097**.
- Patients with more than one expression sample: **118**.
- Distribution of expression samples per patient: 1 sample(s): 979 patients; 2 sample(s): 115 patients; 3 sample(s): 3 patients.
- Expression sample IDs that did not match the documented TCGA barcode pattern: **0**.
- Examples of patients with multiple expression samples (no sample retained or discarded):
  - `TCGA-A7-A0CE`: `TCGA-A7-A0CE-01`, `TCGA-A7-A0CE-11`
  - `TCGA-A7-A0CH`: `TCGA-A7-A0CH-11`, `TCGA-A7-A0CH-01`
  - `TCGA-A7-A0D9`: `TCGA-A7-A0D9-01`, `TCGA-A7-A0D9-11`
  - `TCGA-A7-A0DB`: `TCGA-A7-A0DB-11`, `TCGA-A7-A0DB-01`
  - `TCGA-A7-A0DC`: `TCGA-A7-A0DC-01`, `TCGA-A7-A0DC-11`
  - `TCGA-A7-A13E`: `TCGA-A7-A13E-01`, `TCGA-A7-A13E-11`
  - `TCGA-A7-A13F`: `TCGA-A7-A13F-01`, `TCGA-A7-A13F-11`
  - `TCGA-A7-A13G`: `TCGA-A7-A13G-11`, `TCGA-A7-A13G-01`
  - `TCGA-AC-A23H`: `TCGA-AC-A23H-11`, `TCGA-AC-A23H-01`
  - `TCGA-AC-A2FB`: `TCGA-AC-A2FB-01`, `TCGA-AC-A2FB-11`

## Clinical fields and missingness

Exact field names are discovered from the downloaded clinical header; categories are based on column-name patterns and are not a clinical data dictionary. Multiple receptor fields can reflect different source forms, assays, or specimen contexts; this audit does not select a preferred one.

### Sample ID

`sampleID`, `bcr_sample_barcode`

### Patient ID

`_PATIENT`, `bcr_patient_barcode`, `patient_id`

### PAM50 subtype

`Integrated_Clusters_with_PAM50__nature2012`, `PAM50Call_RNAseq`, `PAM50_mRNA_nature2012`

### ER

`ER_Status_nature2012`, `breast_carcinoma_estrogen_receptor_status`, `breast_carcinoma_immunohistochemistry_er_pos_finding_scale`, `er_detection_method_text`, `er_level_cell_percentage_category`, `metastatic_breast_carcinoma_estrogen_receptor_detection_mthd_txt`, `metastatic_breast_carcinoma_estrogen_receptor_status`, `metastatic_breast_carcinoma_immunohistochemistry_er_pos_cell_scr`, `positive_finding_estrogen_receptor_other_measurement_scale_text`

### PR

`PR_Status_nature2012`, `breast_carcinoma_immunohistochemistry_prgstrn_rcptr_ps_fndng_scl`, `breast_carcinoma_progesterone_receptor_status`, `metastatic_breast_carcinm_ps_fndng_prgstrn_rcptr_thr_msr_scl_txt`, `metastatic_breast_carcinoma_immunohistochemistry_pr_pos_cell_scr`, `metastatic_breast_carcinoma_progesterone_receptor_dtctn_mthd_txt`, `metastatic_breast_carcinoma_progesterone_receptor_status`, `mtsttc_brst_crcnm_mmnhstchmstry_prgstrn_rcptr_pstv_fndng_scl_typ`, `pgr_detection_method_text`, `pos_finding_progesterone_receptor_other_measurement_scale_text`, `progesterone_receptor_level_cell_percent_category`

### HER2

`HER2_Final_Status_nature2012`, `her2_and_centromere_17_positive_finding_other_measuremnt_scl_txt`, `her2_erbb_method_calculation_method_text`, `her2_erbb_pos_finding_cell_percent_category`, `her2_erbb_pos_finding_fluorescence_n_st_hybrdztn_clcltn_mthd_txt`, `her2_immunohistochemistry_level_result`, `her2_neu_and_centromere_17_copy_number_analysis_npt_ttl_nmbr_cnt`, `her2_neu_breast_carcinoma_copy_analysis_input_total_number`, `her2_neu_chromosone_17_signal_ratio_value`, `her2_neu_metastatic_breast_carcinoma_copy_analysis_inpt_ttl_nmbr`, `lab_proc_her2_neu_immunohistochemistry_receptor_status`, `lab_procedure_her2_neu_in_situ_hybrid_outcome_type`, `metastatic_breast_carcinoma_erbb2_immunohistochemistry_levl_rslt`, `metastatic_breast_carcinoma_her2_erbb_method_calculatin_mthd_txt`, `metastatic_breast_carcinoma_her2_erbb_pos_findng_cll_prcnt_ctgry`, `metastatic_breast_carcinoma_her2_neu_chromosone_17_signal_rat_vl`, `pos_finding_her2_erbb2_other_measurement_scale_text`

### Age

`Age_at_Initial_Pathologic_Diagnosis_nature2012`, `age_at_initial_pathologic_diagnosis`, `days_to_birth`, `year_of_initial_pathologic_diagnosis`

### Race/ethnicity

No matching column names are present in this file.

### Survival

`Days_to_Date_of_Last_Contact_nature2012`, `Days_to_date_of_Death_nature2012`, `OS_Time_nature2012`, `OS_event_nature2012`, `Survival_Data_Form_nature2012`, `Vital_Status_nature2012`, `days_to_death`, `days_to_last_followup`, `days_to_last_known_alive`, `vital_status`

Missingness treats blank strings and explicit null-style tokens (`blank`, `--`, `n/a`, `na`, `nan`, `none`, `not applicable`, `not available`, `null`) as missing. Values such as `Indeterminate` and `0` are retained as observed values.

| Category | Exact source column | Present | Missing |
|---|---|---:|---:|
| Sample ID | `sampleID` | 1247/1247 (100.0%) | 0/1247 (0.0%) |
| Sample ID | `bcr_sample_barcode` | 1054/1247 (84.5%) | 193/1247 (15.5%) |
| Patient ID | `_PATIENT` | 1247/1247 (100.0%) | 0/1247 (0.0%) |
| Patient ID | `bcr_patient_barcode` | 1242/1247 (99.6%) | 5/1247 (0.4%) |
| Patient ID | `patient_id` | 1242/1247 (99.6%) | 5/1247 (0.4%) |
| PAM50 subtype | `Integrated_Clusters_with_PAM50__nature2012` | 348/1247 (27.9%) | 899/1247 (72.1%) |
| PAM50 subtype | `PAM50Call_RNAseq` | 956/1247 (76.7%) | 291/1247 (23.3%) |
| PAM50 subtype | `PAM50_mRNA_nature2012` | 522/1247 (41.9%) | 725/1247 (58.1%) |
| ER | `ER_Status_nature2012` | 782/1247 (62.7%) | 465/1247 (37.3%) |
| ER | `breast_carcinoma_estrogen_receptor_status` | 1181/1247 (94.7%) | 66/1247 (5.3%) |
| ER | `breast_carcinoma_immunohistochemistry_er_pos_finding_scale` | 157/1247 (12.6%) | 1090/1247 (87.4%) |
| ER | `er_detection_method_text` | 236/1247 (18.9%) | 1011/1247 (81.1%) |
| ER | `er_level_cell_percentage_category` | 504/1247 (40.4%) | 743/1247 (59.6%) |
| ER | `metastatic_breast_carcinoma_estrogen_receptor_detection_mthd_txt` | 4/1247 (0.3%) | 1243/1247 (99.7%) |
| ER | `metastatic_breast_carcinoma_estrogen_receptor_status` | 49/1247 (3.9%) | 1198/1247 (96.1%) |
| ER | `metastatic_breast_carcinoma_immunohistochemistry_er_pos_cell_scr` | 10/1247 (0.8%) | 1237/1247 (99.2%) |
| ER | `positive_finding_estrogen_receptor_other_measurement_scale_text` | 264/1247 (21.2%) | 983/1247 (78.8%) |
| PR | `PR_Status_nature2012` | 781/1247 (62.6%) | 466/1247 (37.4%) |
| PR | `breast_carcinoma_immunohistochemistry_prgstrn_rcptr_ps_fndng_scl` | 151/1247 (12.1%) | 1096/1247 (87.9%) |
| PR | `breast_carcinoma_progesterone_receptor_status` | 1179/1247 (94.5%) | 68/1247 (5.5%) |
| PR | `metastatic_breast_carcinm_ps_fndng_prgstrn_rcptr_thr_msr_scl_txt` | 4/1247 (0.3%) | 1243/1247 (99.7%) |
| PR | `metastatic_breast_carcinoma_immunohistochemistry_pr_pos_cell_scr` | 11/1247 (0.9%) | 1236/1247 (99.1%) |
| PR | `metastatic_breast_carcinoma_progesterone_receptor_dtctn_mthd_txt` | 4/1247 (0.3%) | 1243/1247 (99.7%) |
| PR | `metastatic_breast_carcinoma_progesterone_receptor_status` | 47/1247 (3.8%) | 1200/1247 (96.2%) |
| PR | `mtsttc_brst_crcnm_mmnhstchmstry_prgstrn_rcptr_pstv_fndng_scl_typ` | 15/1247 (1.2%) | 1232/1247 (98.8%) |
| PR | `pgr_detection_method_text` | 229/1247 (18.4%) | 1018/1247 (81.6%) |
| PR | `pos_finding_progesterone_receptor_other_measurement_scale_text` | 243/1247 (19.5%) | 1004/1247 (80.5%) |
| PR | `progesterone_receptor_level_cell_percent_category` | 458/1247 (36.7%) | 789/1247 (63.3%) |
| HER2 | `HER2_Final_Status_nature2012` | 776/1247 (62.2%) | 471/1247 (37.8%) |
| HER2 | `her2_and_centromere_17_positive_finding_other_measuremnt_scl_txt` | 5/1247 (0.4%) | 1242/1247 (99.6%) |
| HER2 | `her2_erbb_method_calculation_method_text` | 90/1247 (7.2%) | 1157/1247 (92.8%) |
| HER2 | `her2_erbb_pos_finding_cell_percent_category` | 228/1247 (18.3%) | 1019/1247 (81.7%) |
| HER2 | `her2_erbb_pos_finding_fluorescence_n_st_hybrdztn_clcltn_mthd_txt` | 55/1247 (4.4%) | 1192/1247 (95.6%) |
| HER2 | `her2_immunohistochemistry_level_result` | 670/1247 (53.7%) | 577/1247 (46.3%) |
| HER2 | `her2_neu_and_centromere_17_copy_number_analysis_npt_ttl_nmbr_cnt` | 112/1247 (9.0%) | 1135/1247 (91.0%) |
| HER2 | `her2_neu_breast_carcinoma_copy_analysis_input_total_number` | 122/1247 (9.8%) | 1125/1247 (90.2%) |
| HER2 | `her2_neu_chromosone_17_signal_ratio_value` | 249/1247 (20.0%) | 998/1247 (80.0%) |
| HER2 | `her2_neu_metastatic_breast_carcinoma_copy_analysis_inpt_ttl_nmbr` | 3/1247 (0.2%) | 1244/1247 (99.8%) |
| HER2 | `lab_proc_her2_neu_immunohistochemistry_receptor_status` | 1035/1247 (83.0%) | 212/1247 (17.0%) |
| HER2 | `lab_procedure_her2_neu_in_situ_hybrid_outcome_type` | 451/1247 (36.2%) | 796/1247 (63.8%) |
| HER2 | `metastatic_breast_carcinoma_erbb2_immunohistochemistry_levl_rslt` | 21/1247 (1.7%) | 1226/1247 (98.3%) |
| HER2 | `metastatic_breast_carcinoma_her2_erbb_method_calculatin_mthd_txt` | 3/1247 (0.2%) | 1244/1247 (99.8%) |
| HER2 | `metastatic_breast_carcinoma_her2_erbb_pos_findng_cll_prcnt_ctgry` | 11/1247 (0.9%) | 1236/1247 (99.1%) |
| HER2 | `metastatic_breast_carcinoma_her2_neu_chromosone_17_signal_rat_vl` | 9/1247 (0.7%) | 1238/1247 (99.3%) |
| HER2 | `pos_finding_her2_erbb2_other_measurement_scale_text` | 16/1247 (1.3%) | 1231/1247 (98.7%) |
| Age | `Age_at_Initial_Pathologic_Diagnosis_nature2012` | 952/1247 (76.3%) | 295/1247 (23.7%) |
| Age | `age_at_initial_pathologic_diagnosis` | 1242/1247 (99.6%) | 5/1247 (0.4%) |
| Age | `days_to_birth` | 1226/1247 (98.3%) | 21/1247 (1.7%) |
| Age | `year_of_initial_pathologic_diagnosis` | 1240/1247 (99.4%) | 7/1247 (0.6%) |
| Survival | `Days_to_Date_of_Last_Contact_nature2012` | 883/1247 (70.8%) | 364/1247 (29.2%) |
| Survival | `Days_to_date_of_Death_nature2012` | 135/1247 (10.8%) | 1112/1247 (89.2%) |
| Survival | `OS_Time_nature2012` | 952/1247 (76.3%) | 295/1247 (23.7%) |
| Survival | `OS_event_nature2012` | 952/1247 (76.3%) | 295/1247 (23.7%) |
| Survival | `Survival_Data_Form_nature2012` | 952/1247 (76.3%) | 295/1247 (23.7%) |
| Survival | `Vital_Status_nature2012` | 952/1247 (76.3%) | 295/1247 (23.7%) |
| Survival | `days_to_death` | 202/1247 (16.2%) | 1045/1247 (83.8%) |
| Survival | `days_to_last_followup` | 1095/1247 (87.8%) | 152/1247 (12.2%) |
| Survival | `days_to_last_known_alive` | 2/1247 (0.2%) | 1245/1247 (99.8%) |
| Survival | `vital_status` | 1242/1247 (99.6%) | 5/1247 (0.4%) |

## PAM50Call_RNAseq

- Target field: `PAM50Call_RNAseq`.
- Nonmissing labels: **956**; missing labels: **291 / 1,247**.
- Raw label counts (values preserved exactly; no normalization applied):
  - `Basal`: 142
  - `Her2`: 67
  - `LumA`: 434
  - `LumB`: 194
  - `Normal`: 119
- No within-field label pairs differing only by case or whitespace were detected.
- A separate source field, `PAM50_mRNA_nature2012`, uses a different call (522 nonmissing; 725 missing): `Basal-like`: 98, `HER2-enriched`: 58, `Luminal A`: 231, `Luminal B`: 127, `Normal-like`: 8. Do not merge it with `PAM50Call_RNAseq` without a documented decision.
- The selected field uses abbreviated source labels (for example, `LumA`, `Her2`, `Basal`, `Normal`). Possible display-name correspondences to review against project terminology are `LumA`/Luminal A, `LumB`/Luminal B, `Her2`/HER2-enriched, `Basal`/Basal-like, and `Normal`/Normal-like. These are not recodes and have not been applied.

## Important observations and limitations

- Expression and clinical sample counts differ; the exact intersection and all unmatched sample IDs are shown above. Matching is identifier-only, not evidence that labels are valid for every expression profile.
- Patient-level grouping uses the first three TCGA barcode fields. The clinical field named `patient_id` is a short suffix where checked; do not treat it as a unique full barcode without confirmation.
- The source contains multiple PAM50-related fields with different completeness and label vocabularies. `PAM50Call_RNAseq` is reported as-is; no label mapping or sample inclusion decision has been made.
- The clinical file has no column whose name indicates race or ethnicity; demographic subgroup analysis needs a separately documented source/decision.
- HiSeqV2 values are `log2(normalized_count + 1)`, not raw integer counts. They may be used for visualization or documented ML preprocessing, but must not be passed directly to DESeq2/PyDESeq2 as raw counts.
- This is an identifier and metadata audit only. No genes or samples were filtered, no duplicate aliquot was selected, and no train/test split, differential expression, GSEA, or machine learning was performed.

## Decisions reserved for human review

1. Which PAM50 field and label vocabulary should define the later endpoint? This audit does not substitute `PAM50_mRNA_nature2012` for `PAM50Call_RNAseq`.
2. Which ER/PR/HER2 fields should be authoritative if receptor-status analysis is later needed? Candidate fields have different missingness and contexts.
3. How should unmatched samples and multiple samples per patient be handled later? No unmatched IDs or duplicate-patient samples were dropped here.
4. A raw-count source is required before count-based differential expression; HiSeqV2 must remain labeled as normalized expression.

## Reproducibility

Regenerate this report from the repository root with `python scripts/audit_data.py`. The script reads the matrix one row at a time, creates no split, and writes this Markdown report; downloaded source files remain under Git-ignored `data/raw/`.
