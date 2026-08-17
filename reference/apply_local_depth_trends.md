# Apply Local Depth Trends

Enhanced version with Module 8 error handling.

## Usage

``` r
apply_local_depth_trends(
  cokey_data,
  local_predictions,
  unique_depths,
  preserve_correlations = TRUE,
  verbose = getOption("ssurgo.verbose", FALSE)
)
```

## Arguments

- cokey_data:

  Simulation data

- local_predictions:

  Local GP predictions

- unique_depths:

  Depth vector

- preserve_correlations:

  Whether to preserve correlations

- verbose:

  Logical; if `TRUE`, temporarily raises the package's log level so
  `INFO`-level progress messages print for the duration of this call
  (default `FALSE` - quiet). See
  [`set_verbose_logging()`](https://jjmaynard.github.io/soilSIM/reference/set_verbose_logging.md).

## Value

Adjusted simulation data
