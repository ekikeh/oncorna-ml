# Model card — OncoRNA-ML

**Status:** not trained. This scaffold contains no fitted model and makes no performance claims.

## Intended use

Educational/research reproduction of a PAM50 subtype-classification workflow using public, de-identified breast-cancer data.

## Out-of-scope use

Diagnosis, prognosis for an individual, treatment selection, clinical triage, or any patient-care decision. This project is not a medical device and has not been clinically validated.

## Planned evaluation

- Target: PAM50 intrinsic subtype; exact label source and class handling are pending.
- Unit of separation: patient, not aliquot/sample.
- Validation: frozen patient-level test set plus nested cross-validation on training data.
- Primary metric: macro-F1. Additional metrics will include class-wise precision/recall, balanced accuracy, confusion matrix, and calibration/uncertainty if supported by the data.
- Baseline: dummy classifier, followed by regularized multinomial logistic regression.

## Data and limitations

The data source, sample inclusion rules, demographic subgroup availability, label provenance, class imbalance, and preprocessing are not finalized. TCGA is a research cohort and is not population-representative. PAM50 labels inferred from expression can share signal with expression predictors; findings will be framed as label reconstruction, not independent biomarker discovery.

## Update policy

Complete this card after the first locked evaluation. Report dataset version, split definition, software version, results with uncertainty, subgroup limitations, known failure modes, and model artifact checksum. Do not publish a clinical performance claim.
