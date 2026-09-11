# Resolve a property's distribution family, refining `dist = "auto"` from the AOI's own percentile skew rather than a fixed per-property label.

Computes a skewness proxy directly from the AOI's own percentile
triplet: `((high-median)-(median-low))/(high-low)`, aggregated once
across the AOI (default `median`) - NOT per-cell, which would create map
discontinuity artifacts at family-boundary cells.

## Usage

``` r
resolve_property_dist(
  property_config,
  prior_value_rasters,
  prior_probs,
  aggregate_fun = stats::median
)
```

## Arguments

- property_config:

  List with `dist` and optionally `bounds`/`auto_skew_threshold`.

- prior_value_rasters, prior_probs:

  Prior percentile-value rasters/probabilities.

- aggregate_fun:

  Function to aggregate the AOI-wide skew proxy (default
  [`stats::median`](https://rdrr.io/r/stats/median.html)).

## Value

`list(dist=, dist_source= one of "config"/"auto", skew_proxy= NA_real_ unless dist_source=="auto")`.

## Details

Deliberately narrow: only ever resolves to "beta" (if `bounds` is
configured), "normal", or "lognormal" - never "metalog" or "gamma",
since distinguishing those from a 3-point proxy isn't reliable. NEVER
overrides an explicit (non-"auto") `dist` value.

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
[`percentiles_from_draws()`](https://jjmaynard.github.io/soilSIM/reference/percentiles_from_draws.md),
[`quantile_metalog_linear_with_fallback_raster()`](https://jjmaynard.github.io/soilSIM/reference/quantile_metalog_linear_with_fallback_raster.md),
[`remarginalize_awc()`](https://jjmaynard.github.io/soilSIM/reference/remarginalize_awc.md),
[`remarginalize_ensemble_to_posterior()`](https://jjmaynard.github.io/soilSIM/reference/remarginalize_ensemble_to_posterior.md),
[`wrap_nested_rasters()`](https://jjmaynard.github.io/soilSIM/reference/wrap_nested_rasters.md)
