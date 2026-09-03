# Analytic Beta/Gamma percentiles, raster-native, via per-cell `terra::app()`

Unlike the Normal case,
[`qbeta()`](https://rdrr.io/r/stats/Beta.html)/[`qgamma()`](https://rdrr.io/r/stats/GammaDist.html)'s
quantile is a genuinely nonlinear function of both the target
probability AND the shape parameters, so there is no scalar-arithmetic
shortcut - this dispatches per cell via
[`terra::app()`](https://rspatial.github.io/terra/reference/app.html),
but still computes every requested probability for a cell in one
[`vapply()`](https://rdrr.io/r/base/lapply.html) call (R's
[`qbeta()`](https://rdrr.io/r/stats/Beta.html)/[`qgamma()`](https://rdrr.io/r/stats/GammaDist.html)
are already vectorized over their shape-parameter arguments), not one
[`terra::app()`](https://rspatial.github.io/terra/reference/app.html)
pass per probability.

## Usage

``` r
closed_form_percentiles_raster(param1, param2, qfun, posterior_probs)
```

## Arguments

- param1, param2:

  The two shape-parameter SpatRasters (`alpha`/`beta` or
  `shape`/`rate`).

- qfun:

  [`stats::qbeta`](https://rdrr.io/r/stats/Beta.html) or
  [`stats::qgamma`](https://rdrr.io/r/stats/GammaDist.html) (or any
  `function(p, param1, param2)`).

- posterior_probs:

  Numeric probabilities (0-1).

## Value

Named list of percentile `SpatRaster`s.
