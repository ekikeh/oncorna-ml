# Reproducibility appendix

## Historical scientific execution

The source and [data provenance](data_provenance.md), [cohort definition](cohort_definition.md), [frozen split](data_splits.md), [preprocessing contract](preprocessing_v1.md), [training CV](../reports/phase_2b_training_cv.md), and [validation report](../reports/phase_2c_validation.md) document the actual sequence. The split was 506 training, 169 validation, and 169 test patients. The model input was the original ordered 20,530-gene schema with existing `log2(normalized_count + 1)` expression. One four-fold training CV selected logistic `C=1`; it was not a nested CV experiment. The final fit used 675 training-plus-validation patients and retained 17,623 genes, distinct from the 506-patient validation fit retaining 17,613.

Historical reports record versions where established, including Python 3.11.9, NumPy 2.4.6, pandas 2.3.3, scikit-learn 1.9.1, and joblib 1.6.0 for the Phase 2C validation run. These are historical observations for that run; they must not be silently asserted for other phases. Source and artifact hashes appear in the relevant phase records. Missing historical versions are not reconstructed by guessing. Fitted joblib binaries should only be loaded from an authenticated, trusted local source in a compatible environment.

## Completed evaluation evidence

The independently verified, one-time final test is closed. The approved aggregate figures are published in the [Phase 3D report](../reports/phase_3d_aggregate_results.md): 169 test patients, 144 correct classifications, macro F1 0.7934, balanced accuracy 0.7556, and accuracy 0.8521. The public report transcribes approved aggregate information; it is not a recipe to repeat, reopen, or rescore the test. Patient-level predictions and sealed evaluation evidence remain outside the public repository.

## Public repository

The public checkout contains code, versioned configurations, source checksum references, documentation, tests, and a 10-sample by 100-gene **source-derived training preview**. `data/demo/` is not synthetic. The [quick start](../README.md) supports package installation and tests without the full cohort. Rebuilding the real cohort and historical data products requires the original public source files and local derived artifacts; the controlled final evaluation also requires the private fitted artifact and authorization process. Public users cannot reproduce its 169-patient result from this checkout alone.

Current `pyproject.toml` dependency ranges and `environment.yml` target a working development environment; they are not an exact lock of the historical scientific runtime. Existing CI runs lint and tests. Its result should not be interpreted as reproducing the frozen scientific evaluation.

## Standalone synthetic demonstration and planned packaging

The separate [synthetic demonstration](synthetic_demo.md) now generates 100 training and 40 evaluation samples over 12 invented genes with default seed 20261026. It fits its own model in memory and prints aggregate results labeled `synthetic_demo_only`, dataset digests, and dependency versions. It does not load any real fitted artifact or access patient data. Strict input validation is independently reusable; it performs no inference or preprocessing fit.

Exact dependency pins and broader synthetic-only CI/distribution work remain Phase 4D. Any Docker test environment and its dependency set describe the tested software fixtures, not the historical scientific environment. Synthetic scores do not stand in for the final held-out statistics. See [release scope](release_scope.md).
