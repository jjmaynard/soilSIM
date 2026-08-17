# Validate Distribution Fidelity

Enhanced distribution validation utilities.

## Usage

``` r
validate_distribution_fidelity(
  simulation_data,
  simulation_metadata,
  criteria = NULL,
  verbose = getOption("ssurgo.verbose", FALSE)
)
```

## Arguments

- simulation_data:

  Simulation results

- simulation_metadata:

  Metadata from Monte Carlo generation

- criteria:

  Distribution validation criteria

- verbose:

  Logical; if `TRUE`, temporarily raises the package's log level so
  `INFO`-level progress messages print for the duration of this call
  (default `FALSE` - quiet). See
  [`set_verbose_logging()`](https://jjmaynard.github.io/soilSIM/reference/set_verbose_logging.md).

## Value

Distribution fidelity assessment
