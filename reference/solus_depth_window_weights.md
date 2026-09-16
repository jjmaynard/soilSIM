# Trapezoidal Depth-Average Weights over Native SOLUS Depth Slices

Returns the weights that turn SOLUS100's point predictions at the native
depth slices into an estimate of the **depth-average** of the property
over `[top_depth, bottom_depth]` - `sum(weights * slice_values)`
approximates `(1 / (b - a)) * integral_a^b f(d) dd`, where `f` is taken
piecewise-linear between native slices (with `rule = 2` constant
extrapolation outside the `[0, 150]` span SOLUS predicts over). This is
the raster analogue of the thickness-weighted window mean the SSURGO
prior side already computes
([`aggregate_depth_window_by_replicate()`](https://jjmaynard.github.io/soilSIM/reference/aggregate_depth_window_by_replicate.md)),
so the two sides of
[`fuse_property_adaptive()`](https://jjmaynard.github.io/soilSIM/reference/fuse_property_adaptive.md)
represent the same estimand.

## Usage

``` r
solus_depth_window_weights(
  top_depth,
  bottom_depth,
  available_slices = c(0, 5, 15, 30, 60, 100, 150)
)
```

## Arguments

- top_depth, bottom_depth:

  Numeric depth window bounds in cm (`bottom_depth > top_depth`).

- available_slices:

  Native SOLUS depth points (default `c(0, 5, 15, 30, 60, 100, 150)`).

## Value

A named numeric vector over `sort(unique(available_slices))` (names are
the slice depths as character), summing to 1. Zero for slices that do
not contribute.

## Details

Exact for a linear depth profile (reduces to the midpoint value,
matching
[`closest_solus_depth_slice()`](https://jjmaynard.github.io/soilSIM/reference/closest_solus_depth_slice.md));
unbiased to second order for curved profiles (OC decay, clay bulge, pH
trend). Cost is a single weighted raster-stack sum -
`length(available_slices)` is always small, so this scales with the
number of slices, not the number of pixels.

Windows lying entirely at/beyond the deepest slice (or above the
shallowest) fall back to the single nearest slice (constant
extrapolation).

## See also

[`closest_solus_depth_slice()`](https://jjmaynard.github.io/soilSIM/reference/closest_solus_depth_slice.md),
[`fetch_solus_low_pred_high()`](https://jjmaynard.github.io/soilSIM/reference/fetch_solus_low_pred_high.md)

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
[`quantile_metalog_linear_with_fallback_raster()`](https://jjmaynard.github.io/soilSIM/reference/quantile_metalog_linear_with_fallback_raster.md),
[`remarginalize_awc()`](https://jjmaynard.github.io/soilSIM/reference/remarginalize_awc.md),
[`remarginalize_ensemble_to_posterior()`](https://jjmaynard.github.io/soilSIM/reference/remarginalize_ensemble_to_posterior.md),
[`resolve_property_dist()`](https://jjmaynard.github.io/soilSIM/reference/resolve_property_dist.md),
[`simulate_ssurgo_mapunit_draws()`](https://jjmaynard.github.io/soilSIM/reference/simulate_ssurgo_mapunit_draws.md),
[`wrap_nested_rasters()`](https://jjmaynard.github.io/soilSIM/reference/wrap_nested_rasters.md)
