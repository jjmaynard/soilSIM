# Store a value in the cache under `key`

Store a value in the cache under `key`

## Usage

``` r
cache_set(key, kind, value)
```

## Arguments

- key:

  A cache key from
  [`build_cache_key()`](https://jjmaynard.github.io/soilSIM/reference/build_cache_key.md).

- kind:

  Unused; documents intent at call sites. The value is stored keyed only
  by `key`, since
  [`build_cache_key()`](https://jjmaynard.github.io/soilSIM/reference/build_cache_key.md)
  already encodes `kind`.

- value:

  The R object to cache (anything
  [`saveRDS()`](https://rdrr.io/r/base/readRDS.html) can serialize,
  including
  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  objects via their `wrap()`-compatible representation - see Known
  limitation below).

## Value

`TRUE` on success, `FALSE` on a write failure (never errors).

## Known limitation

[`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)/`SpatVector`
objects hold an external pointer to in-memory/on-disk GDAL state that
does not survive a plain
[`saveRDS()`](https://rdrr.io/r/base/readRDS.html)/[`readRDS()`](https://rdrr.io/r/base/readRDS.html)
round-trip as-is. Callers caching raster values must
[`terra::wrap()`](https://rspatial.github.io/terra/reference/wrap.html)
them first and
[`terra::unwrap()`](https://rspatial.github.io/terra/reference/wrap.html)
on read
([`cache_get()`](https://jjmaynard.github.io/soilSIM/reference/cache_get.md)/`cache_set()`
do not do this automatically, since not every cached value is a raster -
e.g. Monte Carlo draw data frames are not).

## See also

Other raster-fusion:
[`align_percentile_probs()`](https://jjmaynard.github.io/soilSIM/reference/align_percentile_probs.md),
[`build_cache_key()`](https://jjmaynard.github.io/soilSIM/reference/build_cache_key.md),
[`cache_get()`](https://jjmaynard.github.io/soilSIM/reference/cache_get.md),
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
[`resolve_property_dist()`](https://jjmaynard.github.io/soilSIM/reference/resolve_property_dist.md),
[`wrap_nested_rasters()`](https://jjmaynard.github.io/soilSIM/reference/wrap_nested_rasters.md)
