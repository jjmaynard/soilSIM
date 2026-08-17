# Optimize GP Hyperparameters

Performs real `n_folds`-fold cross-validation comparing a small set of
[`GPfit::GP_fit()`](https://rdrr.io/pkg/GPfit/man/GP_fit.html)
correlation-family specs (exponential power=1.95, Matern nu=3/2, Matern
nu=5/2 - the practically tunable surface GPfit exposes), picks the
lowest mean cross-validated RMSE, and refits on the full data with the
winning spec. CV results (folds requested, candidates compared, winning
candidate, per-candidate mean RMSE) are attached as the `cv_results`
attribute on the returned model so the search performed is inspectable
without changing the model's class/shape.

## Usage

``` r
optimize_gp_hyperparameters(X, Y, n_folds = 5, gp_control = c(20, 10, 2))
```

## Arguments

- X:

  Scaled predictor matrix (depths).

- Y:

  Response vector.

- n_folds:

  Number of cross-validation folds.

- gp_control:

  Passed through to every
  [`GPfit::GP_fit()`](https://rdrr.io/pkg/GPfit/man/GP_fit.html) call
  this function makes
  ([`k_fold_gp_cv()`](https://jjmaynard.github.io/soilSIM/reference/k_fold_gp_cv.md)'s
  per-fold fits, the winning-candidate refit, and the fallback/baseline
  fits) - see
  [`fit_individual_gp_model()`](https://jjmaynard.github.io/soilSIM/reference/fit_individual_gp_model.md)'s
  docs for why the default is much smaller than `GP_fit()`'s own
  default. With the default 5 folds x 3 correlation candidates + a
  refit, this is ~16 `GP_fit()` calls per invocation - the single
  highest-multiplier fix in the package's performance audit (see
  PERFORMANCE_IMPROVEMENT_PLAN.md).

## Value

A `"GP"`-classed
[`GPfit::GP_fit()`](https://rdrr.io/pkg/GPfit/man/GP_fit.html) model,
with a `cv_results` attribute describing the cross-validation that
selected it.
