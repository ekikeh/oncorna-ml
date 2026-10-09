# Phase 2A — machine-learning experiment design

**Status:** design for human review, 2026-10-09. No classifier has been trained and no model performance has been measured. This document does not authorize opening the frozen validation or test partitions for evaluation.

## 1. Prediction task and scientific scope

Reconstruct the five existing TCGA-BRCA `PAM50Call_RNAseq` labels from the selected primary-tumor RNA expression sample for each of 844 approved patients. The five targets are **Luminal A, Luminal B, Basal-like, HER2-enriched, and Normal-like**, using the fixed normalization in [the cohort definition](cohort_definition.md). This is single-label, five-class classification; every eligible patient has exactly one target and one selected sample. Predictions should include one class per patient and, where the estimator supports them, five aligned class probabilities or decision scores. Store class order explicitly. Probabilities from logistic regression are *model scores*, not calibrated clinical risks; no calibration or decision threshold is part of this first experiment.

The input is the unfiltered `preprocessing_v1` samples-by-genes view: **20,530** source gene rows, retained in source order and identified by the exact source gene-ID strings in `all_gene_schema.json`. Values are already `log2(normalized_count + 1)`; there is no second log transform. `patient_id` and `sample_id` align expression to the frozen split and label table but never become predictors. No clinical fields, subtype metadata, QC flags, or split identifiers enter the feature matrix. The historical 17,623-gene Phase 1E matrix and the full-training 17,613-gene list are provenance artifacts, **not CV inputs**, because their masks were fitted outside individual CV training folds.

The five-class distribution is uneven: the frozen training set has 253 Luminal A, 115 Luminal B, 85 Basal-like, 40 HER2-enriched, and only **13 Normal-like** patients. Accuracy alone would hide rare-class failures. Patient is the independent allocation unit; the current cohort has one sample per patient. If additional aliquots are ever introduced, group them by patient before any resampling.

PAM50 calls are themselves derived from RNA expression. A model that reconstructs these calls may exploit some of the same signal used to assign them. Even perfect agreement would not constitute independent biomarker discovery, an independent biological ground truth, prospective prediction, or clinical validation. This is a reproducible label-reconstruction exercise on a selected TCGA cohort.

## 2. Initial model sequence

| Stage | Model | Purpose and constraints |
|---|---|---|
| 0 | Majority-class dummy classifier | Training-fold class-frequency reference. Fit inside each fold; a deterministic most-frequent prediction is sufficient. It supplies no biological interpretation. |
| 1 | L2-regularized multinomial logistic regression | **Recommended first real classifier.** Joint five-class scores, coefficients tied to source gene IDs, and shrinkage suited to many correlated predictors. Coefficients are descriptive associations with reconstructed labels, not discovered biomarkers. |
| 2, only if justified after Stage 1 review | Linear support vector classifier | A limited linear-margin comparator if Stage 1 exposes a concrete failure or research question. Its decision scores are not probabilities without separate, fold-contained calibration; adding it expands the selection search. |

There are roughly 20,530 candidate genes for 506 training patients before filtering. Correlated genes and the small Normal-like class make variance and selection bias central concerns. Start with the dummy baseline and one regularized linear model. Random forests, boosting, neural networks, and deep learning add complexity and tuning choices without a demonstrated need in this first experiment. Do not add them automatically.

## 3. Leakage-safe development and tuning protocol

1. Load and verify the frozen `split_v1` and `preprocessing_v1` manifests and hashes. Select **only the 506 training patient/sample pairs** for development. Check row-label alignment, the five class counts, unique patient IDs, full ordered gene schema, and finite values before fitting. The split seed `20261006` is historical provenance; use a separately recorded CV seed, proposed `20261009`.
2. Use one **4-fold stratified patient-level CV** partition of the 506 training patients, with shuffling and fixed seed. Each held-out fold has about 3–4 Normal-like patients; each fitting fold has about 9–10. Four folds balance rare-class representation against training size. Persist fold patient IDs and verify no overlap. Do not use the full-training selected-gene list to define the CV input.
3. For every fold and candidate, fit a fresh pipeline on that fold's fitting patients only: `GenePrevalenceFilter(expression_gt=1, minimum_sample_fraction=0.20)` → any imputer if missing values appear → feature scaling → classifier. The filter uses `ceil(0.20 × n_fit)` and applies its learned mask unchanged to the held-out fold. The local all-gene view currently has zero missing values, so omit the imputer initially; if a later input has missing values, use a training-fold-fitted imputer in this position and version the method. Fit a standard scaler after filtering, on fitting patients only, and apply it to the held-out fold. Do not fit any additional feature selection or PCA in this first design. If introduced later, fit it within the same pipeline and every CV fitting fold.
4. Fit the dummy with the same training-fold boundaries. It does not need gene scaling for its predictions, but keeping the partition and input boundary identical makes comparison clear. Score every held-out fold in the fixed five-class order.
5. Tune only logistic regression's L2 strength, with `C ∈ {0.001, 0.01, 0.1, 1}` and an explicit intercept, compatible multinomial solver, fixed convergence tolerance, and sufficiently high iteration cap. Keep class weighting at `None` for the primary run. Optimize **mean fold macro F1**, breaking exact ties toward smaller `C` (stronger regularization). Record convergence warnings; increase the iteration cap or revisit scaling as a documented technical correction, not a silent model search. Do not tune the prevalence threshold, fold count, scaling scheme, class weights, model family, or class mapping against the same CV scores in this phase.
6. Report the selected configuration and its out-of-fold development results as *selection-biased estimates*. Also show each fold's score and support, mean, standard deviation, and range. These four folds are correlated through overlapping training data; their spread is descriptive variability, not a confidence interval or an independent sample of four studies. Keep candidate-level results so the selection is auditable.
7. After choices are locked, fit the chosen full pipeline on all **506 training patients**. Evaluate the frozen **169-patient validation set once** as a development checkpoint and report its class-wise support. If it prompts a revision, log the reason and changed protocol explicitly, treating validation as used for selection. Keep the **169-patient test set locked** until model family, preprocessing, hyperparameters, analysis plan, and report format are frozen. A later final test assessment should be run once. Decide before that assessment whether the final estimator remains trained on 506 patients or is refitted on the combined 675 train-plus-validation patients with the *already fixed* configuration. If refitted, the validation result describes the earlier 506-trained model and cannot independently evaluate the 675-trained model. Never alter the frozen partition membership.

The 4-fold search needs **16 logistic pipeline fits**, plus four dummy fits; the selected full-training fit adds one. The CV estimate is optimistic because it also chooses `C`. Nested CV would offer a less selection-biased internal performance estimate, but with only 13 Normal-like training patients it would make inner-fold minority support even thinner and multiply the fits. It is not necessary for this limited first-stage comparison given the separately frozen validation and final test sets. Repeatedly trying new grids or models against these sets would erode that protection; log every attempted experiment and stop expansion unless a concrete, prespecified reason is documented.

## 4. Evaluation and reporting

**Primary selection metric:** unweighted macro F1 across all five labels, with a fixed label order and `zero_division=0`. It gives each subtype equal weight, though it can change sharply when a rare patient is misclassified. The dummy baseline is the minimum reference, not a performance target.

Supplement with balanced accuracy (mean class recall), overall accuracy only as context, per-class precision/recall/F1 and support, and a confusion matrix with raw counts and clearly labeled row-normalized percentages. For CV, preserve fold-level results and consider an out-of-fold confusion matrix from one prediction per training patient. Distinguish pooled out-of-fold metrics from the arithmetic mean of fold metrics. For validation and eventual test reporting, include denominators and avoid implying precise rare-class performance: Normal-like has only five patients in each held-out partition, so a single error changes its recall by 20 percentage points. No cross-validation fold spread should be presented as a confidence interval for new cohorts. If uncertainty intervals are later reported, state their calculation, sampling unit, and limits for small classes.

Do not compute metrics in Phase 2A. Do not infer biological importance solely from large coefficients: correlated genes, scaling, label circularity, and fold instability complicate interpretation.

## 5. Implementation architecture and reproducible artifacts

Proposed Phase 2B layout, subject to approval:

| Location | Intended content | Git policy |
|---|---|---|
| `configs/modeling_v1.yaml` | Input manifest paths and expected hashes; five-class order; CV seed/folds; fixed filter/scaler; estimator/solver settings and four `C` values; metric definitions; output version | Commit |
| `src/oncorna/modeling.py`, `scripts/run_training_cv.py`, `scripts/evaluate_model.py` | Fold-safe pipeline construction, schema/split assertions, tuning and later evaluation entry points | Commit with focused tests |
| `data/processed/modeling_v1/` | Fold assignments, per-fold and candidate scores, out-of-fold patient predictions, fitted model, hashes, and machine-readable run manifest | Local only; ignored |
| `reports/generated/modeling_v1/` | Generated figures and detailed patient-level diagnostics | Local only; ignored |
| `reports/phase_2b_*.md`, `docs/model_card.md` | Reviewed aggregate results, limitations, methods, provenance references, and model-card update after proper evaluation | Commit only after review; no patient-level tables |

The run manifest should record Git revision and dirty status, config and source/split/preprocessing hashes, ordered sample and gene-schema hashes, fold membership hashes, class order/support, fitted per-fold gene counts, selected `C`, solver convergence status, seeds, command, timestamp, Python/NumPy/pandas/scikit-learn versions, and hashes of generated result/model files. Save fitted models in a format appropriate for the installed library version; loading a Python model artifact requires a trusted source and matching environment. The existing full data, split, derived matrix, predictions, and model binary remain local and outside Git. Avoid committing patient identifiers or source-derived full-cohort values. The small tracked demo is for schema inspection only, never CV or performance claims.

## 6. Compute and limitations

The existing all-gene float64 array is about **132 MiB** (`844 × 20,530 × 8` bytes). A 506-row training view is about 79 MiB before DataFrame, filtered copies, scaling, solver workspace, and saved results. Expect several hundred MiB of active memory, potentially more than 1 GiB depending on copies and library implementation. Begin with one worker to avoid multiplying memory and benchmark one fold before quoting elapsed time; 16 linear-model fits may take minutes to hours depending on hardware and convergence. This is a capacity estimate, not a measured runtime.

Main risks and constraints:

- **Leakage:** the historical cohort-wide 17,623-gene matrix and full-training 17,613-gene mask cannot precede CV filtering; all learned transforms and tuning belong inside training folds. The frozen validation and test patients never train a development pipeline.
- **Expression-derived target:** reconstructed PAM50 calls share the molecular measurement domain with predictors. No independent biomarker or clinical claim follows.
- **High dimension and imbalance:** about 20,530 candidates and only 13 training Normal-like cases increase overfitting and make class-wise results unstable. Regularization and a small search reduce, but do not eliminate, this risk.
- **Generalizability:** one legacy TCGA cohort, selected primary tumors with available labels, two retained low-expression outliers, and no external cohort cannot establish transportability. No patient-level clinical use is supported.
- **Upstream processing:** the Xena matrix is already normalized; whether normalization used information across eventual split boundaries cannot be independently reconstructed here. The fold-safe pipeline controls only downstream learning. The `>1` in 20% prevalence rule was inherited from prior full-cohort QC, so its historical choice may have been influenced by cohort-wide summaries even though each future mask is fitted within a fold. Disclose this fixed-rule history and do not tune the threshold using validation or test data.

## 7. Decisions for human approval before Phase 2B

1. Approve the **dummy → L2 multinomial logistic regression** sequence, with linear SVC reserved for a separately justified follow-up.
2. Approve one seeded **4-fold stratified CV**, the fixed prevalence filter, training-fold scaling, no imputer for the present zero-missing matrix, the four-value `C` grid, and macro F1 selection with smaller-`C` tie break.
3. Approve the validation checkpoint policy and choose, before final testing, whether the final estimator will train on 506 or refit on 675 patients. The proposed default is to refit on 675 only after the 506-trained validation assessment and all choices are frozen.
4. Confirm that the Normal-like class remains included despite its small support and that the existing primary-tumor selection and retained QC outliers remain unchanged.

This report stops at design. Implementation, fitting, validation evaluation, and test evaluation require a subsequent phase.
