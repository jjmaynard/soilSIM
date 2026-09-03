# Depth-Trend GP Adjustment, Guarded by `GPfit` Availability

Applies
[`apply_local_gp_adjustments()`](https://jjmaynard.github.io/soilSIM/reference/apply_local_gp_adjustments.md)
(`R/multivariate-adjustment.R` - fits its own local GP per property from
each cokey's own within-simulation depth trend, no pre-supplied GP
models needed) per cokey, when `GPfit` is installed and a cokey has
enough distinct depths. Cokeys with fewer than `min_depths` distinct
depths pass through unadjusted, exactly as the original per-cokey guard
did.

## Usage

``` r
maybe_adjust_soil_data_depth_trend(
  sim_long,
  properties,
  min_depths = 2,
  parallel = FALSE,
  n_cores = NULL,
  config = NULL
)
```

## Arguments

- sim_long:

  Long-format simulated data with `cokey`, `hzdept_r`, and property
  columns.

- properties:

  Character vector of property column names to adjust.

- min_depths:

  Minimum distinct depths required to attempt GP fitting (default 2,
  matching the source's `length(unique_depths) >= 2` guard).

- parallel:

  Logical; if `TRUE`, process cokeys across multiple
  [`future::multisession`](https://future.futureverse.org/reference/multisession.html)
  worker processes (default `FALSE` - sequential, matching prior
  behavior exactly). Falls back to sequential processing if the parallel
  setup itself errors. See
  [`run_parallel_lapply()`](https://jjmaynard.github.io/soilSIM/reference/run_parallel_lapply.md)
  (`R/parallel-utils.R`).

- n_cores:

  Number of worker processes to use when `parallel = TRUE` (default
  `max(1, parallel::detectCores() - 1)`).

- config:

  Optional Monte Carlo config, passed through to
  [`adjust_one_cokey_depth_trend()`](https://jjmaynard.github.io/soilSIM/reference/adjust_one_cokey_depth_trend.md)
  -\>
  [`apply_local_gp_adjustments()`](https://jjmaynard.github.io/soilSIM/reference/apply_local_gp_adjustments.md)
  -\>
  [`apply_gp_depth_trends()`](https://jjmaynard.github.io/soilSIM/reference/apply_gp_depth_trends.md).
  `config$monte_carlo$vertical_correlation_method` (default
  `"joint_copula"` as of `VERTICAL_CORRELATION_IMPROVEMENT_PLAN.md`
  Phase 13; set to `"gp_quantile_retrofit"` to opt back into the
  original algorithm) selects the vertical-correlation method; `NULL`
  (default) resolves to `"joint_copula"`, matching
  [`get_monte_carlo_defaults()`](https://jjmaynard.github.io/soilSIM/reference/get_monte_carlo_defaults.md)'s
  own default.

## Value

`sim_long`, depth-trend-adjusted where possible.

## Details

Each cokey's GP fitting is completely independent of every other
cokey's, so this step is embarrassingly parallel - profiling on a real
AOI showed it as the dominant cost of the whole SSURGO simulation
pipeline (see
[`apply_local_gp_adjustments()`](https://jjmaynard.github.io/soilSIM/reference/apply_local_gp_adjustments.md)/[`fit_local_gp_model_single()`](https://jjmaynard.github.io/soilSIM/reference/fit_local_gp_model_single.md)),
so for AOIs with many cokeys, `parallel = TRUE` can give a further
speedup roughly proportional to available cores on top of the
sequential-path optimizations already applied there.
