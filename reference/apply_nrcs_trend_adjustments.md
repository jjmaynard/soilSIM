# Apply NRCS Trend Adjustments

Enhanced version with Module 8 integration and proper Module 5 function
calls.

## Usage

``` r
apply_nrcs_trend_adjustments(
  cokey_data,
  gp_models,
  model_group,
  properties,
  preserve_correlations = TRUE,
  verbose = getOption("ssurgo.verbose", FALSE)
)
```

## Arguments

- cokey_data:

  Simulation data for a single cokey

- gp_models:

  NRCS GP models from gp_modeling module

- model_group:

  GP model group for this cokey

- properties:

  Properties to adjust

- preserve_correlations:

  Whether to preserve correlations

- verbose:

  Logical; if `TRUE`, temporarily raises the package's log level so
  `INFO`-level progress messages print for the duration of this call
  (default `FALSE` - quiet). See
  [`set_verbose_logging()`](https://jjmaynard.github.io/soilSIM/reference/set_verbose_logging.md).

## Value

Adjusted simulation data
