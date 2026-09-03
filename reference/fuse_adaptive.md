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
  [`bayesian_update()`](https://jjmaynard.github.io/soilSIM/reference/bayesian_update.md)
  route is used; above it, the closed-form route.

- n_samples:

  Samples drawn per side for the general route's KDE fusion.

- grid_resolution:

  Passed to
  [`bayesian_update()`](https://jjmaynard.github.io/soilSIM/reference/bayesian_update.md)
  for the general route. `NULL` (default) resolves to
  [FUSE_GENERAL_KDE_DEFAULT_GRID_RESOLUTION](https://jjmaynard.github.io/soilSIM/reference/FUSE_GENERAL_KDE_DEFAULT_GRID_RESOLUTION.md)
  (`0.1`) - see that constant's docs for the profiling/accuracy
  justification for why this differs from
  [`bayesian_update()`](https://jjmaynard.github.io/soilSIM/reference/bayesian_update.md)'s
  own standalone default of `0.01`.

- verbose:

  If TRUE (default), print the chosen route and why.

- mukey_raster, mukey_draws:

  Optional - opts into `prior_fusion_method = "raw_draws"` on whichever
  route runs
  ([`fuse_general_kde()`](https://jjmaynard.github.io/soilSIM/reference/fuse_general_kde.md)
  for the general route - task P2.2/P2.3;
  [`fuse_closed_form()`](https://jjmaynard.github.io/soilSIM/reference/fuse_closed_form.md)
  for the closed-form route - task P2.6, added at negligible cost since
  that route's raw_draws fit is a per-mukey precompute +
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
  [`bayesian_update()`](https://jjmaynard.github.io/soilSIM/reference/bayesian_update.md)'s
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
  `posterior_probs` entry, on both routes (as of P3.3).

- `route`: one of `"bayesian_update_general"`, `"closed_form_normal"`,
  `"closed_form_beta"`, `"closed_form_gamma"`.

- `route_detail`: for `family %in% c("beta","gamma")` on the closed-form
  route, a `SpatRaster` of 1 (native same-family fusion) / 0 (per-cell
  fallback). `NULL` for "normal" and for the general route.

- `n_fallback_cells`: count of cells that used the fallback.

- `diagnostics`:
  `list(ncell=, threshold_cells=, family=, elapsed_sec=)`.
