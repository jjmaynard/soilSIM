# Fetch a Raster of SSURGO Map Unit Keys for an AOI

Fetch a Raster of SSURGO Map Unit Keys for an AOI

## Usage

``` r
fetch_ssurgo_mukey_raster(aoi_vect)
```

## Arguments

- aoi_vect:

  A
  [`terra::SpatVector`](https://rspatial.github.io/terra/reference/SpatVector-class.html)
  (projected, e.g. EPSG:5070) - the AOI footprint itself, not widened to
  its bounding box (unlike `R/adapter-ssurgo-acquire.R`'s
  [`process_aoi_and_get_mukeys_working()`](https://jjmaynard.github.io/soilSIM/reference/process_aoi_and_get_mukeys_working.md),
  which bbox-widens for its tabular by-mukey-list query -
  [`soilDB::mukey.wcs()`](http://ncss-tech.github.io/soilDB/reference/mukey.wcs.md)
  bbox-clips internally regardless of input polygon shape, so both
  resolve to the same extent; passing the footprint here avoids widening
  the *cached* grid's extent unnecessarily for callers that only need
  the AOI itself covered).

## Value

A single-layer, categorical (factor)
[`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
named `"mukey"`, or `NULL` if the grid fetch returns nothing.

## Implementation note

[`soilDB::mukey.wcs()`](http://ncss-tech.github.io/soilDB/reference/mukey.wcs.md)
returns a raster whose cell values are the mukey codes
(`unique(values(res))`), so no separate polygon fetch/rasterize step is
needed - the function issues exactly one network call per cache miss.
The result is disk-cached (`cache.R`, kind `"mukey_grid"`,
depth-agnostic key via
[`mukey_grid_cache_key()`](https://jjmaynard.github.io/soilSIM/reference/mukey_grid_cache_key.md))
so repeated calls for the same AOI across separate top-level user calls
hit zero network calls.

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
[`fetch_ssurgo_percentiles()`](https://jjmaynard.github.io/soilSIM/reference/fetch_ssurgo_percentiles.md),
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
