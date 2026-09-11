# Raw Per-Mukey Monte Carlo Draws, Keyed by Mukey (Not Collapsed to Percentiles)

The raw-draws analogue of
[`percentiles_from_draws`](https://jjmaynard.github.io/soilSIM/reference/percentiles_from_draws.md):
groups `draws` by `mukey` the same way, but keeps each mukey's full
vector of simulated values instead of collapsing it to a handful of
quantiles. This is what lets raster-fusion routes (`R/core-fusion.R`'s
[`fuse_general_kde()`](https://jjmaynard.github.io/soilSIM/reference/fuse_general_kde.md)/[`fuse_texture_group_batch_core()`](https://jjmaynard.github.io/soilSIM/reference/fuse_texture_group_batch.md),
once they opt into `prior_fusion_method = "raw_draws"`) fuse directly
against the real empirical distribution computed once per unique mukey,
rather than resampling a lower-fidelity reconstruction from just a few
percentile values.

## Usage

``` r
lookup_mukey_draws(draws, property_id)
```

## Arguments

- draws:

  A data frame from
  [`simulate_ssurgo_mapunit_draws`](https://jjmaynard.github.io/soilSIM/reference/simulate_ssurgo_mapunit_draws.md),
  for the same AOI/depth window the caller's mukey raster covers.

- property_id:

  One of
  [`property_to_sim_column`](https://jjmaynard.github.io/soilSIM/reference/property_to_sim_column.md)'s
  recognized ids.

## Value

A named list keyed by mukey (as character), each element a numeric
vector of that mukey's simulated values (across every cokey/replicate) -
or `NULL` if `property_id`'s simulated column isn't present in `draws`.
