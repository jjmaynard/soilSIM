# Vectorized feasibility check: a valid metalog quantile function must be monotonically increasing in y (equivalent to its density staying non-negative everywhere). Probes `quantile_metalog_linear_raster()` at a fixed grid of y-values and flags cells where consecutive probe values decrease - a probe, not a proof, but fully vectorized raster arithmetic. Streams one probe raster at a time rather than materializing all of them (avoids an allocation failure at large cell counts).

Vectorized feasibility check: a valid metalog quantile function must be
monotonically increasing in y (equivalent to its density staying
non-negative everywhere). Probes
[`quantile_metalog_linear_raster()`](https://jjmaynard.github.io/soilSIM/reference/fit_metalog_linear_raster.md)
at a fixed grid of y-values and flags cells where consecutive probe
values decrease - a probe, not a proof, but fully vectorized raster
arithmetic. Streams one probe raster at a time rather than materializing
all of them (avoids an allocation failure at large cell counts).

## Usage

``` r
check_metalog_feasibility_raster(
  fit,
  bounds,
  boundedness,
  y_grid = seq(0.02, 0.98, by = 0.02)
)
```

## Arguments

- fit:

  Output of
  [`fit_metalog_linear_raster()`](https://jjmaynard.github.io/soilSIM/reference/fit_metalog_linear_raster.md).

- bounds, boundedness:

  Same as
  [`fit_metalog_linear_raster()`](https://jjmaynard.github.io/soilSIM/reference/fit_metalog_linear_raster.md)'s
  arguments for this `fit`.

- y_grid:

  Probability grid to probe.

## See also

Other raster-fusion:
[`align_percentile_probs()`](https://jjmaynard.github.io/soilSIM/reference/align_percentile_probs.md),
[`build_cache_key()`](https://jjmaynard.github.io/soilSIM/reference/build_cache_key.md),
[`cache_get()`](https://jjmaynard.github.io/soilSIM/reference/cache_get.md),
[`cache_set()`](https://jjmaynard.github.io/soilSIM/reference/cache_set.md),
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
