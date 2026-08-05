# Perform GP Cross-Validation

Runs real k-fold cross-validation (reusing the shared
[`k_fold_gp_cv()`](https://jjmaynard.github.io/soilSIM/reference/k_fold_gp_cv.md)
helper from `gp-modeling.R`) for each property in `gp_models` against
the matching column of `training_data`, aggregating mean CV RMSE and an
R-squared derived from it.

## Usage

``` r
perform_gp_cross_validation(gp_models, training_data, criteria)
```

## Arguments

- gp_models:

  GP models from
  [`build_stratified_gp_models()`](https://jjmaynard.github.io/soilSIM/reference/build_stratified_gp_models.md).

- training_data:

  Original training data with `hzdept_r` and one column per property in
  `gp_models`.

- criteria:

  List, optionally with `n_folds` (default 5).

## Value

List with `cv_rmse`, `cv_r_squared`, `cv_quality`, `cv_method`.
