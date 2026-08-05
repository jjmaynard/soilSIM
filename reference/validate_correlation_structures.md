# Validate Correlation Structures

Enhanced correlation validation utilities and Module 6 integration.

## Usage

``` r
validate_correlation_structures(
  correlation_matrices,
  simulation_data,
  config = NULL,
  verbose = getOption("ssurgo.verbose", FALSE)
)
```

## Arguments

- correlation_matrices:

  Correlation matrices from correlation_structure module

- simulation_data:

  Final simulation data

- config:

  Correlation validation configuration

- verbose:

  Logical; if `TRUE`, temporarily raises the package's log level so
  `INFO`-level progress messages print for the duration of this call
  (default `FALSE` - quiet). See
  [`set_verbose_logging()`](https://jjmaynard.github.io/soilSIM/reference/set_verbose_logging.md).

## Value

Correlation validation results
