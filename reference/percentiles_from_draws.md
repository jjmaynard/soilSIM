# Compute Per-Mukey Percentile Rasters from Already-Simulated SSURGO Draws

The quantile/rasterize half of
[`fetch_ssurgo_percentiles()`](https://jjmaynard.github.io/soilSIM/reference/fetch_ssurgo_percentiles.md),
factored out so a single shared
[`simulate_ssurgo_mapunit_draws()`](https://jjmaynard.github.io/soilSIM/reference/simulate_ssurgo_mapunit_draws.md)
result can be reused across multiple properties instead of resimulating
once per property -
[`simulate_cokey_generalized()`](https://jjmaynard.github.io/soilSIM/reference/simulate_cokey_generalized.md)
already simulates every recognized property jointly in one pass per
cokey, so a caller that needs several properties from the same AOI/depth
window (e.g.
[`run_fusion_group()`](https://jjmaynard.github.io/soilSIM/reference/run_fusion_group.md)'s
texture members) only needs to run the (expensive) simulation once and
call this per property afterward.

## Usage

``` r
percentiles_from_draws(
  mukey_raster,
  draws,
  property_id,
  probs = c(0.05, 0.25, 0.5, 0.75, 0.95)
)
```

## Arguments

- mukey_raster:

  A categorical (factor) mukey
  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  from
  [`fetch_ssurgo_mukey_raster()`](https://jjmaynard.github.io/soilSIM/reference/fetch_ssurgo_mukey_raster.md).

- draws:

  A data frame from
  [`simulate_ssurgo_mapunit_draws()`](https://jjmaynard.github.io/soilSIM/reference/simulate_ssurgo_mapunit_draws.md),
  for the same AOI/depth window `mukey_raster` covers.

- property_id:

  One of
  [`property_to_sim_column()`](https://jjmaynard.github.io/soilSIM/reference/property_to_sim_column.md)'s
  recognized ids.

- probs:

  Percentile probabilities to compute (default
  `c(0.05, 0.25, 0.5, 0.75, 0.95)`).

## Value

`list(values = <named list of percentile-value SpatRasters>, probs = probs)`,
or `NULL` if `property_id`'s simulated column isn't present in `draws`.

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
[`fuse_adaptive()`](https://jjmaynard.github.io/soilSIM/reference/fuse_adaptive.md),
[`quantile_metalog_linear_with_fallback_raster()`](https://jjmaynard.github.io/soilSIM/reference/quantile_metalog_linear_with_fallback_raster.md),
[`remarginalize_awc()`](https://jjmaynard.github.io/soilSIM/reference/remarginalize_awc.md),
[`remarginalize_ensemble_to_posterior()`](https://jjmaynard.github.io/soilSIM/reference/remarginalize_ensemble_to_posterior.md),
[`resolve_property_dist()`](https://jjmaynard.github.io/soilSIM/reference/resolve_property_dist.md),
[`wrap_nested_rasters()`](https://jjmaynard.github.io/soilSIM/reference/wrap_nested_rasters.md)
