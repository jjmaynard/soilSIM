# Fetch SSURGO Percentile-Value Rasters for an AOI

The top-level SSURGO "prior" entry point for `R/core-fusion.R`'s
[`fuse_property_adaptive()`](https://jjmaynard.github.io/soilSIM/reference/fuse_property_adaptive.md):
rasterizes map units, Monte Carlo-simulates the requested property, and
returns per-cell percentile-value rasters in the `list(values=, probs=)`
shape
[`fuse_property_adaptive()`](https://jjmaynard.github.io/soilSIM/reference/fuse_property_adaptive.md)
expects.

## Usage

``` r
fetch_ssurgo_percentiles(
  aoi_vect,
  property_id,
  top_depth,
  bottom_depth,
  probs = c(0.05, 0.25, 0.5, 0.75, 0.95),
  n_mc = 1000,
  parallel = FALSE,
  n_cores = NULL,
  requested_properties = NULL,
  seed = NULL
)
```

## Arguments

- aoi_vect:

  A
  [`terra::SpatVector`](https://rspatial.github.io/terra/reference/SpatVector-class.html)
  AOI.

- property_id:

  One of
  [`property_to_sim_column()`](https://jjmaynard.github.io/soilSIM/reference/property_to_sim_column.md)'s
  recognized ids.

- top_depth, bottom_depth:

  Numeric depth window bounds in cm.

- probs:

  Percentile probabilities to compute (default
  `c(0.05, 0.25, 0.5, 0.75, 0.95)`).

- n_mc:

  Number of triangular draws per component (default 1000).

- parallel, n_cores:

  Passed through to
  [`simulate_ssurgo_mapunit_draws()`](https://jjmaynard.github.io/soilSIM/reference/simulate_ssurgo_mapunit_draws.md)'s
  `parallel`/`n_cores` - the per-cokey depth-trend GP fitting step is
  this function's dominant cost for AOIs with many cokeys. Default
  `parallel = FALSE` matches prior behavior exactly.

- requested_properties:

  Passed through to
  [`simulate_ssurgo_mapunit_draws()`](https://jjmaynard.github.io/soilSIM/reference/simulate_ssurgo_mapunit_draws.md)
  to restrict the simulation. Default `NULL` here means "just
  `property_id`" (this function only ever reads one property's column
  back out), so a bare call simulates only what it needs. Pass an
  explicit vector to keep extra properties in the shared draws (e.g. a
  caller reusing them); pass `character(0)`-safe values only - `NULL` is
  the "just this property" shortcut, not "all".

- seed:

  Optional integer for opt-in determinism, forwarded to
  [`simulate_ssurgo_mapunit_draws()`](https://jjmaynard.github.io/soilSIM/reference/simulate_ssurgo_mapunit_draws.md).
  `NULL` (default) = current stochastic behavior.

## Value

`list(values = <named list of percentile-value SpatRasters>, probs = probs)`,
or `NULL` if the mukey raster or the Monte Carlo draws are unavailable
for this AOI.

## See also

Other raster-fusion:
[`align_percentile_probs()`](https://jjmaynard.github.io/soilSIM/reference/align_percentile_probs.md),
[`build_cache_key()`](https://jjmaynard.github.io/soilSIM/reference/build_cache_key.md),
[`cache_get()`](https://jjmaynard.github.io/soilSIM/reference/cache_get.md),
[`cache_set()`](https://jjmaynard.github.io/soilSIM/reference/cache_set.md),
[`check_metalog_feasibility_raster()`](https://jjmaynard.github.io/soilSIM/reference/check_metalog_feasibility_raster.md),
[`closest_solus_depth_slice()`](https://jjmaynard.github.io/soilSIM/reference/closest_solus_depth_slice.md),
[`fetch_solus_low_pred_high()`](https://jjmaynard.github.io/soilSIM/reference/fetch_solus_low_pred_high.md),
[`fetch_solus_percentiles()`](https://jjmaynard.github.io/soilSIM/reference/fetch_solus_percentiles.md),
[`fetch_solus_percentiles_multiproperty()`](https://jjmaynard.github.io/soilSIM/reference/fetch_solus_percentiles_multiproperty.md),
[`fetch_ssurgo_mukey_raster()`](https://jjmaynard.github.io/soilSIM/reference/fetch_ssurgo_mukey_raster.md),
[`fit_beta_mle_newton_raster()`](https://jjmaynard.github.io/soilSIM/reference/fit_beta_mle_newton_raster.md),
[`fit_beta_mom_raster()`](https://jjmaynard.github.io/soilSIM/reference/fit_beta_mom_raster.md),
[`fit_gamma_mom_raster()`](https://jjmaynard.github.io/soilSIM/reference/fit_gamma_mom_raster.md),
[`fit_metalog_linear_raster()`](https://jjmaynard.github.io/soilSIM/reference/fit_metalog_linear_raster.md),
[`fit_normal_raster()`](https://jjmaynard.github.io/soilSIM/reference/fit_normal_raster.md),
[`fuse_adaptive()`](https://jjmaynard.github.io/soilSIM/reference/fuse_adaptive.md),
[`percentiles_from_draws()`](https://jjmaynard.github.io/soilSIM/reference/percentiles_from_draws.md),
[`property_to_sim_column()`](https://jjmaynard.github.io/soilSIM/reference/property_to_sim_column.md),
[`quantile_metalog_linear_with_fallback_raster()`](https://jjmaynard.github.io/soilSIM/reference/quantile_metalog_linear_with_fallback_raster.md),
[`remarginalize_awc()`](https://jjmaynard.github.io/soilSIM/reference/remarginalize_awc.md),
[`remarginalize_ensemble_to_posterior()`](https://jjmaynard.github.io/soilSIM/reference/remarginalize_ensemble_to_posterior.md),
[`resolve_property_dist()`](https://jjmaynard.github.io/soilSIM/reference/resolve_property_dist.md),
[`simulate_ssurgo_mapunit_draws()`](https://jjmaynard.github.io/soilSIM/reference/simulate_ssurgo_mapunit_draws.md),
[`solus_depth_window_weights()`](https://jjmaynard.github.io/soilSIM/reference/solus_depth_window_weights.md),
[`wrap_nested_rasters()`](https://jjmaynard.github.io/soilSIM/reference/wrap_nested_rasters.md)
