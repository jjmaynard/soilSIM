# Apply GP Depth Trends with Correlation Preservation

Enhanced core function that applies GP-derived depth trends using Module
8 utilities for robust error handling and validation.

## Usage

``` r
apply_gp_depth_trends(
  cokey_data,
  gp_predictions,
  properties,
  preserve_correlations = TRUE,
  primary_property = NULL,
  verbose = getOption("ssurgo.verbose", FALSE)
)
```

## Arguments

- cokey_data:

  Simulation data for a single cokey

- gp_predictions:

  Named list of GP predictions by property

- properties:

  Properties to adjust

- preserve_correlations:

  Whether to preserve correlations

- primary_property:

  Reference property for correlation preservation

- verbose:

  Logical; if `TRUE`, temporarily raises the package's log level so
  `INFO`-level progress messages print for the duration of this call
  (default `FALSE` - quiet). See
  [`set_verbose_logging()`](https://jjmaynard.github.io/soilSIM/reference/set_verbose_logging.md).

## Value

Adjusted simulation data
