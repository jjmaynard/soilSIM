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
optimize_gp_hyperparameters(X, Y, n_folds = 5)
```

## Arguments

- X:

  Scaled predictor matrix (depths).

- Y:

  Response vector.

- n_folds:

  Number of cross-validation folds.

## Value

A `"GP"`-classed
[`GPfit::GP_fit()`](https://rdrr.io/pkg/GPfit/man/GP_fit.html) model,
with a `cv_results` attribute describing the cross-validation that
selected it.
