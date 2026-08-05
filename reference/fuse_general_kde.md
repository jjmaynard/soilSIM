# Fully general fusion path for `fuse_adaptive()`: per cell, draw samples from both sides' percentiles (via `simulate_from_percentiles(method= "linear_cdf")`), fuse via `bayesian_update()`, and moment-match the posterior samples back into the requested family so the output contract matches the closed-form route.

Fully general fusion path for
[`fuse_adaptive()`](https://jjmaynard.github.io/soilSIM/reference/fuse_adaptive.md):
per cell, draw samples from both sides' percentiles (via
`simulate_from_percentiles(method= "linear_cdf")`), fuse via
[`bayesian_update()`](https://jjmaynard.github.io/soilSIM/reference/bayesian_update.md),
and moment-match the posterior samples back into the requested family so
the output contract matches the closed-form route.

## Usage

``` r
fuse_general_kde(
  prior_value_rasters,
  lik_value_rasters,
  percentile_probs,
  family,
  bounds,
  n_samples,
  grid_resolution
)
```

## Arguments

- prior_value_rasters, lik_value_rasters:

  Lists of percentile-value SpatRasters.

- percentile_probs:

  Numeric probabilities matching the value-raster list order.

- family:

  One of "normal", "beta", "gamma".

- bounds:

  Required for `family = "beta"`.

- n_samples:

  Samples drawn per side, per cell.

- grid_resolution:

  Passed to
  [`bayesian_update()`](https://jjmaynard.github.io/soilSIM/reference/bayesian_update.md).
