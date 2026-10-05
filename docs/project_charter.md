# Project charter — OncoRNA-ML

**Status:** draft for Phase 0; finalize data-source details before downloading or analyzing a cohort.

## Research question

Can a reproducible machine-learning workflow reconstruct TCGA-BRCA PAM50 intrinsic-subtype labels from gene-expression measurements, and do model errors or important features align with known breast-cancer biology?

This is a computational learning exercise, not clinical validation and not a search for a deployable diagnostic.

## Planned endpoints

- **Primary:** multiclass PAM50 subtype classification. The label source and rules for rare/unknown classes must be documented before modeling.
- **Secondary (optional):** overall-survival risk analysis, only if time, event, censoring, and patient identifiers can be harmonized and independently checked.
- **Descriptive:** Basal-like versus Luminal A differential expression and pathway enrichment, with a clear statement about which samples were used and why.

## Inclusion and data rules

- Public, de-identified, open-access data only; no controlled BAM/FASTQ and no attempts to re-identify participants.
- Prefer primary tumors for subtype classification; define normal samples and duplicate aliquots before analysis.
- One patient contributes to only one split. Resolve multiple samples per patient using a documented rule before splitting.
- Keep raw counts for negative-binomial differential-expression methods. Never treat TPM/FPKM as raw counts.
- Check redistribution terms before placing any data in Git. Full-cohort and derived matrices stay out of the public repository.

## Validation plan

- Generate and freeze a patient-level test split before examining test performance.
- Keep all preprocessing that learns from data—including scaling, feature selection, and tuning—inside training folds.
- Compare against a dummy baseline and a regularized multinomial logistic-regression model before adding complex models.
- Headline macro-F1; also report balanced accuracy, per-class precision/recall, confusion matrix, and uncertainty where sample size permits.
- Use a synthetic demo only to check that the code runs. Synthetic results are not evidence of biological or clinical performance.

## Success criteria

1. A new user can create the environment and run tests from a clean clone.
2. A deterministic synthetic/demo workflow can later run without downloading the full cohort.
3. Input provenance, patient-level split logic, parameters, and software versions are recorded.
4. Any eventual results include leakage checks, baselines, class-wise errors, limitations, and a clear non-clinical disclaimer.

## Open decisions before Phase 1

- Exact TCGA-BRCA expression source and release, and whether it supplies raw counts or normalized expression.
- PAM50 label source and sample-ID mapping.
- Patient/sample inclusion rules and treatment of duplicate aliquots or rare classes.
- Demo-data generation method and redistribution/license statement.
