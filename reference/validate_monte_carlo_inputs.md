# Validate Monte Carlo Inputs (Enhanced)

Comprehensive input validation

## Usage

``` r
validate_monte_carlo_inputs(
  soil_data,
  properties,
  correlation_matrix,
  n_realizations,
  config,
  verbose = getOption("ssurgo.verbose", FALSE)
)
```

## Arguments

- soil_data:

  Input soil data

- properties:

  Properties to simulate

- correlation_matrix:

  Optional correlation matrix

- n_realizations:

  Number of realizations

- config:

  Configuration settings

- verbose:

  Logical; if `TRUE`, temporarily raises the package's log level so
  `INFO`-level progress messages print for the duration of this call
  (default `FALSE` - quiet). See
  [`set_verbose_logging()`](https://jjmaynard.github.io/soilSIM/reference/set_verbose_logging.md).

## Value

Validation results
