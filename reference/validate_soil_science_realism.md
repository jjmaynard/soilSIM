# Validate Soil Science Realism

Enhanced soil science validation property validation and utilities.

## Usage

``` r
validate_soil_science_realism(
  simulation_data,
  original_data = NULL,
  config = NULL,
  verbose = getOption("ssurgo.verbose", FALSE)
)
```

## Arguments

- simulation_data:

  Final simulation data

- original_data:

  Original SSURGO data for comparison

- config:

  Soil science validation configuration

- verbose:

  Logical; if `TRUE`, temporarily raises the package's log level so
  `INFO`-level progress messages print for the duration of this call
  (default `FALSE` - quiet). See
  [`set_verbose_logging()`](https://jjmaynard.github.io/soilSIM/reference/set_verbose_logging.md).

## Value

Soil science validation results
