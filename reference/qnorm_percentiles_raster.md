# Analytic Normal percentiles, raster-native (`mu + sigma * qnorm(p)`)

Pure raster arithmetic - `qnorm(p)` is a single scalar per requested
probability (unlike
[`qbeta()`](https://rdrr.io/r/stats/Beta.html)/[`qgamma()`](https://rdrr.io/r/stats/GammaDist.html),
whose quantile depends on the shape parameters too, not just `p`), so no
per-cell dispatch
([`terra::app()`](https://rspatial.github.io/terra/reference/app.html))
is needed at all; exact, and effectively free.

## Usage

``` r
qnorm_percentiles_raster(mu, sigma, posterior_probs)
```

## Arguments

- mu, sigma:

  Normal parameter SpatRasters.

- posterior_probs:

  Numeric probabilities (0-1).

## Value

Named list of percentile `SpatRaster`s (names `"P1"`, `"P50"`, ...
matching `posterior_probs`).
