# Extract NRCS Depth Trends

Enhanced version with Module 8 validation and error handling.

## Usage

``` r
extract_nrcs_depth_trends(
  gp_models,
  properties,
  depths = seq(0, 200, by = 10),
  verbose = getOption("ssurgo.verbose", FALSE)
)
```

## Arguments

- gp_models:

  NRCS GP models

- properties:

  Properties to extract trends for

- depths:

  Depths for trend extraction

- verbose:

  Logical; if `TRUE`, temporarily raises the package's log level so
  `INFO`-level progress messages print for the duration of this call
  (default `FALSE` - quiet). See
  [`set_verbose_logging()`](https://jjmaynard.github.io/soilSIM/reference/set_verbose_logging.md).

## Value

List of depth trends by property and group
