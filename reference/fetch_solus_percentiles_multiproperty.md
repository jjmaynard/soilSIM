# Fetch SOLUS100 Percentile-Value Rasters for Multiple Variables, One Depth Window (S1)

Batched sibling of
[`fetch_solus_percentiles()`](https://jjmaynard.github.io/soilSIM/reference/fetch_solus_percentiles.md):
thin per-variable `list(values=, probs=)` packaging over
[`fetch_solus_low_pred_high_multiproperty()`](https://jjmaynard.github.io/soilSIM/reference/fetch_solus_low_pred_high_multiproperty.md)'s
one-call fetch.

## Usage

``` r
fetch_solus_percentiles_multiproperty(
  aoi_vect,
  solus_variables,
  top_depth,
  bottom_depth
)
```

## Arguments

- aoi_vect:

  A
  [`terra::SpatVector`](https://rspatial.github.io/terra/reference/SpatVector-class.html)
  AOI.

- solus_variables:

  Character vector of
  [`soilDB::fetchSOLUS()`](http://ncss-tech.github.io/soilDB/reference/fetchSOLUS.md)-recognized
  variable names.

- top_depth, bottom_depth:

  Numeric depth window bounds in cm.

## Value

A named list, one entry per element of `solus_variables`, each either
`list(values = list(P025=, P50=, P975=), probs = c(0.025, 0.5, 0.975))`
(matching
[`fetch_solus_percentiles()`](https://jjmaynard.github.io/soilSIM/reference/fetch_solus_percentiles.md)'s
shape) or `NULL` if any of that variable's low/pred/high rasters is
unavailable.

## See also

[`fetch_solus_percentiles()`](https://jjmaynard.github.io/soilSIM/reference/fetch_solus_percentiles.md),
[`fetch_solus_low_pred_high_multiproperty()`](https://jjmaynard.github.io/soilSIM/reference/fetch_solus_low_pred_high_multiproperty.md)

Other raster-fusion:
[`align_percentile_probs()`](https://jjmaynard.github.io/soilSIM/reference/align_percentile_probs.md),
[`build_cache_key()`](https://jjmaynard.github.io/soilSIM/reference/build_cache_key.md),
[`cache_get()`](https://jjmaynard.github.io/soilSIM/reference/cache_get.md),
[`cache_set()`](https://jjmaynard.github.io/soilSIM/reference/cache_set.md),
[`check_metalog_feasibility_raster()`](https://jjmaynard.github.io/soilSIM/reference/check_metalog_feasibility_raster.md),
[`closest_solus_depth_slice()`](https://jjmaynard.github.io/soilSIM/reference/closest_solus_depth_slice.md),
[`fetch_solus_low_pred_high()`](https://jjmaynard.github.io/soilSIM/reference/fetch_solus_low_pred_high.md),
[`fetch_solus_percentiles()`](https://jjmaynard.github.io/soilSIM/reference/fetch_solus_percentiles.md),
[`fit_beta_mle_newton_raster()`](https://jjmaynard.github.io/soilSIM/reference/fit_beta_mle_newton_raster.md),
[`fit_beta_mom_raster()`](https://jjmaynard.github.io/soilSIM/reference/fit_beta_mom_raster.md),
[`fit_gamma_mom_raster()`](https://jjmaynard.github.io/soilSIM/reference/fit_gamma_mom_raster.md),
[`fit_metalog_linear_raster()`](https://jjmaynard.github.io/soilSIM/reference/fit_metalog_linear_raster.md),
[`fit_normal_raster()`](https://jjmaynard.github.io/soilSIM/reference/fit_normal_raster.md),
[`fuse_adaptive()`](https://jjmaynard.github.io/soilSIM/reference/fuse_adaptive.md),
[`percentiles_from_draws()`](https://jjmaynard.github.io/soilSIM/reference/percentiles_from_draws.md),
[`quantile_metalog_linear_with_fallback_raster()`](https://jjmaynard.github.io/soilSIM/reference/quantile_metalog_linear_with_fallback_raster.md),
[`remarginalize_awc()`](https://jjmaynard.github.io/soilSIM/reference/remarginalize_awc.md),
[`remarginalize_ensemble_to_posterior()`](https://jjmaynard.github.io/soilSIM/reference/remarginalize_ensemble_to_posterior.md),
[`resolve_property_dist()`](https://jjmaynard.github.io/soilSIM/reference/resolve_property_dist.md),
[`solus_depth_window_weights()`](https://jjmaynard.github.io/soilSIM/reference/solus_depth_window_weights.md),
[`wrap_nested_rasters()`](https://jjmaynard.github.io/soilSIM/reference/wrap_nested_rasters.md)
