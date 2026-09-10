# Closed-form same-family fusion path for `fuse_adaptive()`, with a per-cell fallback to Normal-moment fusion (re-expressed back into the requested family) where the naive same-family fusion is infeasible.

Closed-form same-family fusion path for
[`fuse_adaptive()`](https://jjmaynard.github.io/soilSIM/reference/fuse_adaptive.md),
with a per-cell fallback to Normal-moment fusion (re-expressed back into
the requested family) where the naive same-family fusion is infeasible.

## Usage

``` r
fuse_closed_form(
  prior_value_rasters,
  lik_value_rasters,
  percentile_probs,
  family,
  bounds,
  posterior_probs = NULL,
  mukey_raster = NULL,
  mukey_draws = NULL
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

- posterior_probs:

  Numeric probabilities (0-1) at which to report posterior percentiles.
  `NULL` (default) resolves to
  [FUSE_POSTERIOR_DEFAULT_PROBS](https://jjmaynard.github.io/soilSIM/reference/FUSE_POSTERIOR_DEFAULT_PROBS.md).
  Computed analytically (`qnorm`/`qbeta`/`qgamma` at the fused family
  parameters) - exact regardless of how extreme the requested
  probability is, unlike
  [`fuse_general_kde()`](https://jjmaynard.github.io/soilSIM/reference/fuse_general_kde.md)'s
  sample/grid-based route.

- mukey_raster, mukey_draws:

  Optional - opts into `prior_fusion_method = "raw_draws"` for this
  (closed-form) route: fits the PRIOR side's family parameters directly
  from each cell's mukey's real Monte Carlo draws (via
  [`mukey_draws_closed_form_fit_raster()`](https://jjmaynard.github.io/soilSIM/reference/mukey_draws_closed_form_fit_raster.md))
  instead of
  [`fit_normal_raster()`](https://jjmaynard.github.io/soilSIM/reference/fit_normal_raster.md)/
  [`fit_beta_mle_newton_raster()`](https://jjmaynard.github.io/soilSIM/reference/fit_beta_mle_newton_raster.md)/[`fit_gamma_mom_raster()`](https://jjmaynard.github.io/soilSIM/reference/fit_gamma_mom_raster.md)'s
  percentile-triplet formulas - cells whose mukey has no draws entry
  fall back to the percentile-based fit for that cell only, same
  convention as
  [`fuse_general_kde()`](https://jjmaynard.github.io/soilSIM/reference/fuse_general_kde.md)'s
  raw_draws branch. The LIKELIHOOD (SOLUS) side is unaffected either
  way - no raw-draws equivalent exists for it. `NULL` (default) for
  either parameter preserves the original percentile-only behavior
  exactly.

## Value

`list(posterior=, route_detail=, n_fallback_cells=)`. `posterior`
carries a new `percentiles` field (named list of `SpatRaster`s)
alongside the existing family-native parameters.
