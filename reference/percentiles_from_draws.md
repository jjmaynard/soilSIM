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
[`run_stage1_fusion_group()`](https://jjmaynard.github.io/soilSIM/reference/run_stage1_fusion_group.md)'s
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
