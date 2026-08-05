# Apply Local GP Adjustments

Enhanced version with Module 8 utilities and better error handling.

## Usage

``` r
apply_local_gp_adjustments(
  cokey_data,
  properties,
  preserve_correlations = TRUE,
  min_depths = 3,
  config = NULL,
  verbose = getOption("ssurgo.verbose", FALSE)
)
```

## Arguments

- cokey_data:

  Simulation data for a single cokey

- properties:

  Properties to adjust

- preserve_correlations:

  Whether to preserve correlations

- min_depths:

  Minimum depths required for GP fitting

- config:

  Configuration from Module 8

- verbose:

  Logical; if `TRUE`, temporarily raises the package's log level so
  `INFO`-level progress messages print for the duration of this call
  (default `FALSE` - quiet). See
  [`set_verbose_logging()`](https://jjmaynard.github.io/soilSIM/reference/set_verbose_logging.md).

## Value

Adjusted simulation data
