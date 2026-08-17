# Fit a Beta distribution's (alpha, beta) rasters via vectorized Newton-Raphson MLE.

If a percentile value lands exactly on `bounds` (rescaled x = 0 or 1
exactly), `log(x)`/`log(1-x)` in the score equations would be `-Inf`,
degenerating the fit to the clamp floor for every cell uniformly -
unlikely with real empirical percentiles but a real risk for synthetic
or explicitly bound-recorded data. Clamped away from the exact boundary
by a tiny epsilon - a no-op for realistic input.

## Usage

``` r
fit_beta_mle_newton_raster(value_rasters, bounds, n_iter = 15, eps = 1e-06)

quantile_beta_mle_newton_raster(fit, q)
```

## Arguments

- value_rasters:

  List of percentile-value SpatRasters (any count \>= 3).

- bounds:

  c(lower, upper) physical bounds.

- n_iter:

  Fixed Newton-Raphson iteration count (validated sufficient at 15
  upstream).

- eps:

  Clamp epsilon away from the exact 0/1 boundary.

- fit:

  Output of `fit_beta_mle_newton_raster()`.

- q:

  Target quantile probability.

## Value

list(alpha = SpatRaster, beta = SpatRaster)
