# Project charter — OncoRNA-ML

**Status:** Data-preparation Phases 1A–1G are implemented. The selected public source, approved cohort, processed matrix, frozen patient-level split, and compact attributed demo are documented in `docs/`; no model training or evaluation, differential expression, enrichment analysis, or biological interpretation has started.

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
- A compact, fixed demo subset from the training partition is committed for table/schema inspection; synthetic fixtures may be used in tests. Neither the demo nor synthetic fixtures are evidence of biological or clinical performance.

## Success criteria

1. A new user can create the environment and run tests from a clean clone.
2. A new user can inspect the committed demo without downloading the full cohort; reproducing it is deterministic when the documented local Phase 1C/1E/1F inputs are present.
3. Input provenance, patient-level split logic, parameters, and software versions are recorded.
4. Any eventual results include leakage checks, baselines, class-wise errors, limitations, and a clear non-clinical disclaimer.

## Recorded Phase 1 decisions and future review

- Cohort/sample inclusion, patient-to-sample mapping, PAM50 normalization, QC, matrix semantics, and the immutable patient-level split are recorded in the corresponding Phase 1 documentation.
- Phase 1G documents and tests the small training-only demo sampling rules and source attribution. Xena does not display a file-level license for the source files; review current upstream terms before any broader redistribution. The repository's MIT license covers project code only.
- If count-based differential expression is implemented later, select a raw-count source; the Phase 1A Xena expression file is log2-normalized and is not suitable as raw DESeq2/PyDESeq2 input.
- Before any Phase 2 modeling or evaluation, follow the locked-test policy in `docs/data_splits.md` and keep learned preprocessing, feature selection, and tuning inside training folds.
