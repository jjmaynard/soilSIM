# Validate Simulation Output (Enhanced)

Enhanced output validation

## Usage

``` r
validate_simulation_output(
  simulation_results,
  properties,
  config,
  verbose = getOption("ssurgo.verbose", FALSE)
)
```

## Arguments

- simulation_results:

  Simulation results array

- properties:

  Property names

- config:

  Configuration settings

- verbose:

  Logical; if `TRUE`, temporarily raises the package's log level so
  `INFO`-level progress messages print for the duration of this call
  (default `FALSE` - quiet). See
  [`set_verbose_logging()`](https://jjmaynard.github.io/soilSIM/reference/set_verbose_logging.md).

## Value

Enhanced validation results
