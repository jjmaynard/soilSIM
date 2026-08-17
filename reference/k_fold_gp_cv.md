# K-Fold Cross-Validation for GP Correlation Structure Selection

Runs `n_folds`-fold cross-validation of
[`GPfit::GP_fit()`](https://rdrr.io/pkg/GPfit/man/GP_fit.html) for each
supplied correlation-family candidate, returning the mean held-out RMSE
per candidate. Shared helper used by
[`optimize_gp_hyperparameters()`](https://jjmaynard.github.io/soilSIM/reference/optimize_gp_hyperparameters.md)
and
[`perform_gp_cross_validation()`](https://jjmaynard.github.io/soilSIM/reference/perform_gp_cross_validation.md)
(validation-diagnostics.R) so the fold-splitting logic is not duplicated
across files.

## Usage

``` r
k_fold_gp_cv(X, Y, n_folds, corr_candidates, gp_control = c(20, 10, 2))
```

## Arguments

- X:

  Predictor matrix (one row per observation).

- Y:

  Response vector.

- n_folds:

  Number of cross-validation folds.

- corr_candidates:

  Named list of `corr` specs to pass to
  [`GPfit::GP_fit()`](https://rdrr.io/pkg/GPfit/man/GP_fit.html) (e.g.
  `list(type = "exponential", power = 1.95)`).

- gp_control:

  Passed through to each fold's
  [`GPfit::GP_fit()`](https://rdrr.io/pkg/GPfit/man/GP_fit.html)
  `control` argument - see
  [`fit_individual_gp_model()`](https://jjmaynard.github.io/soilSIM/reference/fit_individual_gp_model.md)'s
  docs for why the default is much smaller than `GP_fit()`'s own
  default. This is the highest-multiplier call site in the package:
  `n_folds x length(corr_candidates)` `GP_fit()` calls per invocation.

## Value

List with `mean_rmse_by_candidate`, a named list of mean held-out RMSE
per candidate (`NA` for a candidate that failed on every fold).
