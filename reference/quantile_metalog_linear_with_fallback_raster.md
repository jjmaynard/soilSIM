# Metalog quantile with automatic fallback to `linear_cdf` for infeasible cells - zero effect on feasible cells, exact `linear_cdf` match on infeasible ones.

Feasibility is computed ONCE per fit (via
[`check_metalog_feasibility_raster()`](https://jjmaynard.github.io/soilSIM/reference/check_metalog_feasibility_raster.md))
and passed in rather than recomputed per call, since it depends only on
the fitted coefficients, not on which quantile `q` is being evaluated.

## Usage

``` r
quantile_metalog_linear_with_fallback_raster(
  fit,
  infeasible_r,
  full_value_rasters,
  full_probs,
  q,
  bounds,
  boundedness
)
```

## Arguments

- fit:

  Output of
  [`fit_metalog_linear_raster()`](https://jjmaynard.github.io/soilSIM/reference/fit_metalog_linear_raster.md)
  (fit on INTERIOR percentiles).

- infeasible_r:

  Output of
  [`check_metalog_feasibility_raster()`](https://jjmaynard.github.io/soilSIM/reference/check_metalog_feasibility_raster.md)
  for this same `fit`.

- full_value_rasters, full_probs:

  The FULL percentile set (including p=0/p=1 if available) -
  `linear_cdf` has no restriction there and benefits from every
  available knot.

- q:

  Target quantile probability.

- bounds, boundedness:

  Same as
  [`fit_metalog_linear_raster()`](https://jjmaynard.github.io/soilSIM/reference/fit_metalog_linear_raster.md)'s
  arguments.

## Value

list(value = blended quantile raster, used_fallback = infeasible_r
verbatim)

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
[`fetch_ssurgo_percentiles()`](https://jjmaynard.github.io/soilSIM/reference/fetch_ssurgo_percentiles.md),
[`fit_beta_mle_newton_raster()`](https://jjmaynard.github.io/soilSIM/reference/fit_beta_mle_newton_raster.md),
[`fit_beta_mom_raster()`](https://jjmaynard.github.io/soilSIM/reference/fit_beta_mom_raster.md),
[`fit_gamma_mom_raster()`](https://jjmaynard.github.io/soilSIM/reference/fit_gamma_mom_raster.md),
[`fit_metalog_linear_raster()`](https://jjmaynard.github.io/soilSIM/reference/fit_metalog_linear_raster.md),
[`fit_normal_raster()`](https://jjmaynard.github.io/soilSIM/reference/fit_normal_raster.md),
[`fuse_adaptive()`](https://jjmaynard.github.io/soilSIM/reference/fuse_adaptive.md),
[`percentiles_from_draws()`](https://jjmaynard.github.io/soilSIM/reference/percentiles_from_draws.md),
[`property_to_sim_column()`](https://jjmaynard.github.io/soilSIM/reference/property_to_sim_column.md),
[`remarginalize_awc()`](https://jjmaynard.github.io/soilSIM/reference/remarginalize_awc.md),
[`remarginalize_ensemble_to_posterior()`](https://jjmaynard.github.io/soilSIM/reference/remarginalize_ensemble_to_posterior.md),
[`resolve_property_dist()`](https://jjmaynard.github.io/soilSIM/reference/resolve_property_dist.md),
[`simulate_ssurgo_mapunit_draws()`](https://jjmaynard.github.io/soilSIM/reference/simulate_ssurgo_mapunit_draws.md),
[`solus_depth_window_weights()`](https://jjmaynard.github.io/soilSIM/reference/solus_depth_window_weights.md),
[`wrap_nested_rasters()`](https://jjmaynard.github.io/soilSIM/reference/wrap_nested_rasters.md)
