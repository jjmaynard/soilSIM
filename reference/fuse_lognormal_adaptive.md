# Lognormal adapter for `fuse_property_adaptive()`: size-adaptive general/closed-form route on log-transformed values.

Lognormal adapter for
[`fuse_property_adaptive()`](https://jjmaynard.github.io/soilSIM/reference/fuse_property_adaptive.md):
size-adaptive general/closed-form route on log-transformed values.

## Usage

``` r
fuse_lognormal_adaptive(
  prior_value_rasters,
  prior_probs,
  lik_value_rasters,
  lik_probs,
  ncell,
  threshold_cells,
  n_samples = 500,
  grid_resolution = NULL,
  verbose = TRUE,
  posterior_probs = NULL,
  mukey_raster = NULL,
  mukey_draws = NULL
)
```

## Arguments

- posterior_probs:

  Numeric probabilities (0-1) at which to report posterior percentiles.
  `NULL` (default) resolves to
  [FUSE_POSTERIOR_DEFAULT_PROBS](https://jjmaynard.github.io/soilSIM/reference/FUSE_POSTERIOR_DEFAULT_PROBS.md).
  Small-AOI branch: forwarded to
  [`fuse_adaptive()`](https://jjmaynard.github.io/soilSIM/reference/fuse_adaptive.md)'s
  general route unchanged (see that function's docs). Large-AOI
  closed-form branch: computed via
  [`qlnorm()`](https://rdrr.io/r/stats/Lognormal.html) on the LOG-SPACE
  fused parameters (`posterior_log$mu`/`sigma`, captured before the
  moment-match-back-to-raw-space step below) - more faithful to the
  assumed lognormal family than deriving percentiles from the raw-space
  Normal-moment approximation (`posterior$mu`/`sigma`) would be.

- mukey_raster, mukey_draws:

  Optional raw-draws inputs (see
  [`fuse_adaptive()`](https://jjmaynard.github.io/soilSIM/reference/fuse_adaptive.md)'s
  docs). Small-AOI branch: forwarded to
  [`fuse_adaptive()`](https://jjmaynard.github.io/soilSIM/reference/fuse_adaptive.md)'s
  general route unchanged. Large-AOI closed-form branch Raw-draws
  inputs: each mukey's real draws are fit directly in log-space
  (`mukey_draws_closed_form_fit_raster(..., family = "lognormal")`) and
  merged over the percentile-triplet-derived log-space fit wherever a
  mukey has draws coverage - the same merge pattern used for
  normal/beta/gamma on the closed-form route.
