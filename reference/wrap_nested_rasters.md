# Wrap/unwrap every `SpatRaster` found anywhere inside an arbitrarily-nested list

A generic counterpart to
[`wrap_percentile_list()`](https://jjmaynard.github.io/soilSIM/reference/wrap_percentile_list.md)/[`unwrap_percentile_list()`](https://jjmaynard.github.io/soilSIM/reference/wrap_percentile_list.md),
for result structures whose shape isn't the fixed
`list(values=, probs=)` percentile form - e.g.
[`run_fusion()`](https://jjmaynard.github.io/soilSIM/reference/run_fusion.md)'s
return value, which nests `SpatRaster`s at `prior$values`,
`likelihood$values`, and inside `posterior` (whose own shape varies by
`dist`). Recurses into every list element, replacing each `SpatRaster`
with
[`terra::wrap()`](https://rspatial.github.io/terra/reference/wrap.html)'s
`PackedSpatRaster` representation (or the reverse via
[`terra::unwrap()`](https://rspatial.github.io/terra/reference/wrap.html))
and leaving everything else untouched - directly addresses
[`cache_set()`](https://jjmaynard.github.io/soilSIM/reference/cache_set.md)'s
own `@section Known limitation:` for callers caching a full fusion
result rather than a raw percentile fetch.

## Usage

``` r
wrap_nested_rasters(x)

unwrap_nested_rasters(x)
```

## Arguments

- x:

  A list, arbitrarily nested, that may contain `SpatRaster` objects at
  any depth.

## Value

`x`, with every `SpatRaster` replaced by its wrapped/unwrapped
equivalent.

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
[`resolve_property_dist()`](https://jjmaynard.github.io/soilSIM/reference/resolve_property_dist.md)
