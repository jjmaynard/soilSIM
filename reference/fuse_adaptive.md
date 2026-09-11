# Fuse a prior and likelihood belief distribution, choosing the fusion route by AOI size rather than requiring the caller to pick.

Fuse a prior and likelihood belief distribution, choosing the fusion
route by AOI size rather than requiring the caller to pick.

## Usage

``` r
fuse_adaptive(
  prior_value_rasters,
  lik_value_rasters,
  percentile_probs,
  family = c("normal", "beta", "gamma"),
  bounds = NULL,
  threshold_cells = 80000,
  n_samples = 500,
  grid_resolution = NULL,
  verbose = TRUE,
  mukey_raster = NULL,
  mukey_draws = NULL,
  posterior_probs = NULL
)
```

## Arguments

- prior_value_rasters, lik_value_rasters:

  Named lists of percentile-value `SpatRaster`s (e.g.
  `list(P0=, P5=, P50=, P95=, P100=)`), same layer count/order on both
  sides.

- percentile_probs:

  Numeric probabilities matching the value-raster list order (e.g.
  `c(0, 0.05, 0.5, 0.95, 1)`).

- family:

  One of "normal", "beta", "gamma" - the family both sides are fit in
  for the closed-form route, and the family the general route's output
  is moment-matched back into for a consistent output contract across
  routes.

- bounds:

  Required for `family = "beta"` (`c(lower, upper)`).

- threshold_cells:

  AOI cell count at or below which the general
  [`update_prior()`](https://jjmaynard.github.io/soilSIM/reference/update_prior.md)
  route is used; above it, the closed-form route.

- n_samples:

  Samples drawn per side for the general route's KDE fusion.

- grid_resolution:

  Passed to
  [`update_prior()`](https://jjmaynard.github.io/soilSIM/reference/update_prior.md)
  for the general route. `NULL` (default) resolves to
  [FUSE_GENERAL_KDE_DEFAULT_GRID_RESOLUTION](https://jjmaynard.github.io/soilSIM/reference/FUSE_GENERAL_KDE_DEFAULT_GRID_RESOLUTION.md)
  (`0.1`) - see that constant's docs for the profiling/accuracy
  justification for why this differs from
  [`update_prior()`](https://jjmaynard.github.io/soilSIM/reference/update_prior.md)'s
  own standalone default of `0.01`.

- verbose:

  If TRUE (default), print the chosen route and why.

- mukey_raster, mukey_draws:

  Optional - opts into `prior_fusion_method = "raw_draws"` on whichever
  route runs
  ([`fuse_general_kde()`](https://jjmaynard.github.io/soilSIM/reference/fuse_general_kde.md)
  for the general route -;
  [`fuse_closed_form()`](https://jjmaynard.github.io/soilSIM/reference/fuse_closed_form.md)
  for the closed-form route -, added at negligible cost since that
  route's raw_draws fit is a per-mukey precompute +
  [`terra::subst()`](https://rspatial.github.io/terra/reference/subst.html)
  broadcast, not a per-cell computation - see
  [`mukey_draws_closed_form_fit_raster()`](https://jjmaynard.github.io/soilSIM/reference/mukey_draws_closed_form_fit_raster.md)'s
  docs). `NULL` (default) preserves original behavior exactly on both
  routes.

- posterior_probs:

  Passed through to whichever route runs
  ([`fuse_general_kde()`](https://jjmaynard.github.io/soilSIM/reference/fuse_general_kde.md)
  for the general route,
  [`fuse_closed_form()`](https://jjmaynard.github.io/soilSIM/reference/fuse_closed_form.md)
  for the closed-form route - see
  [FUSE_POSTERIOR_DEFAULT_PROBS](https://jjmaynard.github.io/soilSIM/reference/FUSE_POSTERIOR_DEFAULT_PROBS.md)).
  `NULL` (default) resolves to the same rich default on both routes. The
  general route computes percentiles from
  [`update_prior()`](https://jjmaynard.github.io/soilSIM/reference/update_prior.md)'s
  discretized posterior (exact relative to `grid_resolution`, but tail
  percentiles are limited by how well the KDE approximates the tails
  from a finite sample); the closed-form route computes them
  analytically (`qnorm`/`qbeta`/ `qgamma` at the fused family
  parameters - exact regardless of how extreme the probability).

## Value

A list:

- `posterior`: family-native parameter list (`mu`/`sigma` for "normal",
  `alpha`/`beta` for "beta", `shape`/`rate` for "gamma") - always this
  shape regardless of which route ran. Also carries `percentiles`: a
  named list of `SpatRaster`s (`"P01"`..`"P99"` by default), one per
  `posterior_probs` entry, on both routes.

- `route`: one of `"bayesian_update_general"`, `"closed_form_normal"`,
  `"closed_form_beta"`, `"closed_form_gamma"`.

- `route_detail`: for `family %in% c("beta","gamma")` on the closed-form
  route, a `SpatRaster` of 1 (native same-family fusion) / 0 (per-cell
  fallback). `NULL` for "normal" and for the general route.

- `n_fallback_cells`: count of cells that used the fallback.

- `diagnostics`:
  `list(ncell=, threshold_cells=, family=, elapsed_sec=)`.

## See also

Other raster-fusion:
[`align_percentile_probs()`](https://jjmaynard.github.io/soilSIM/reference/align_percentile_probs.md),
[`build_cache_key()`](https://jjmaynard.github.io/soilSIM/reference/build_cache_key.md),
[`cache_get()`](https://jjmaynard.github.io/soilSIM/reference/cache_get.md),
[`cache_set()`](https://jjmaynard.github.io/soilSIM/reference/cache_set.md),
[`check_metalog_feasibility_raster()`](https://jjmaynard.github.io/soilSIM/reference/check_metalog_feasibility_raster.md),
[`fit_beta_mle_newton_raster()`](https://jjmaynard.github.io/soilSIM/reference/fit_beta_mle_newton_raster.md),
[`fit_beta_mom_raster()`](https://jjmaynard.github.io/soilSIM/reference/fit_beta_mom_raster.md),
[`fit_gamma_mom_raster()`](https://jjmaynard.github.io/soilSIM/reference/fit_gamma_mom_raster.md),
[`fit_metalog_linear_raster()`](https://jjmaynard.github.io/soilSIM/reference/fit_metalog_linear_raster.md),
[`fit_normal_raster()`](https://jjmaynard.github.io/soilSIM/reference/fit_normal_raster.md),
[`percentiles_from_draws()`](https://jjmaynard.github.io/soilSIM/reference/percentiles_from_draws.md),
[`quantile_metalog_linear_with_fallback_raster()`](https://jjmaynard.github.io/soilSIM/reference/quantile_metalog_linear_with_fallback_raster.md),
[`remarginalize_awc()`](https://jjmaynard.github.io/soilSIM/reference/remarginalize_awc.md),
[`remarginalize_ensemble_to_posterior()`](https://jjmaynard.github.io/soilSIM/reference/remarginalize_ensemble_to_posterior.md),
[`resolve_property_dist()`](https://jjmaynard.github.io/soilSIM/reference/resolve_property_dist.md),
[`wrap_nested_rasters()`](https://jjmaynard.github.io/soilSIM/reference/wrap_nested_rasters.md)
