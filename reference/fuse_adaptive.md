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
  verbose = TRUE
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
  for the general route.

- verbose:

  If TRUE (default), print the chosen route and why.

## Value

A list:

- `posterior`: family-native parameter list (`mu`/`sigma` for "normal",
  `alpha`/`beta` for "beta", `shape`/`rate` for "gamma") - always this
  shape regardless of which route ran.

- `route`: one of `"bayesian_update_general"`, `"closed_form_normal"`,
  `"closed_form_beta"`, `"closed_form_gamma"`.

- `route_detail`: for `family %in% c("beta","gamma")` on the closed-form
  route, a `SpatRaster` of 1 (native same-family fusion) / 0 (per-cell
  fallback). `NULL` for "normal" and for the general route.

- `n_fallback_cells`: count of cells that used the fallback.

- `diagnostics`:
  `list(ncell=, threshold_cells=, family=, elapsed_sec=)`.
