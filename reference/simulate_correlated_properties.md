# Simulate Correlated Properties Using Enhanced Distribution Framework

Core simulation function leveraging Module 0 statistical utilities

## Usage

``` r
simulate_correlated_properties(
  simulation_params,
  correlation_matrix,
  n_realizations,
  properties,
  config,
  verbose = getOption("ssurgo.verbose", FALSE)
)
```

## Arguments

- simulation_params:

  List containing simulation parameters for each horizon

- correlation_matrix:

  Correlation matrix for the properties

- n_realizations:

  Number of realizations to generate

- properties:

  Character vector of property names

- config:

  Configuration settings

- verbose:

  Logical; if `TRUE`, temporarily raises the package's log level so
  `INFO`-level progress messages print for the duration of this call
  (default `FALSE` - quiet). See
  [`set_verbose_logging()`](https://jjmaynard.github.io/soilSIM/reference/set_verbose_logging.md).

## Value

Array of simulated values, dimensioned horizons by properties by
realizations
