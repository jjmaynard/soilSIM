# Assess Property Constraints

Enhanced property constraint validation property utilities.

## Usage

``` r
assess_property_constraints(
  simulation_data,
  criteria = NULL,
  verbose = getOption("ssurgo.verbose", FALSE)
)
```

## Arguments

- simulation_data:

  Simulation data

- criteria:

  Property constraint criteria

- verbose:

  Logical; if `TRUE`, temporarily raises the package's log level so
  `INFO`-level progress messages print for the duration of this call
  (default `FALSE` - quiet). See
  [`set_verbose_logging()`](https://jjmaynard.github.io/soilSIM/reference/set_verbose_logging.md).

## Value

Property constraint assessment
