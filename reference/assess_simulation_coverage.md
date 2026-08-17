# Assess Simulation Coverage

Enhanced coverage assessment validation utilities.

## Usage

``` r
assess_simulation_coverage(
  simulation_data,
  original_data,
  criteria = NULL,
  verbose = getOption("ssurgo.verbose", FALSE)
)
```

## Arguments

- simulation_data:

  Simulation results

- original_data:

  Original data for comparison

- criteria:

  Coverage assessment criteria

- verbose:

  Logical; if `TRUE`, temporarily raises the package's log level so
  `INFO`-level progress messages print for the duration of this call
  (default `FALSE` - quiet). See
  [`set_verbose_logging()`](https://jjmaynard.github.io/soilSIM/reference/set_verbose_logging.md).

## Value

Coverage assessment results
