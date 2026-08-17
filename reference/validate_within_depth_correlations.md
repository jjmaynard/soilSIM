# Validate Within Depth Correlations

Enhanced within-depth correlation validation utilities.

## Usage

``` r
validate_within_depth_correlations(
  simulation_data,
  criteria = NULL,
  verbose = getOption("ssurgo.verbose", FALSE)
)
```

## Arguments

- simulation_data:

  Simulation data

- criteria:

  Depth validation criteria

- verbose:

  Logical; if `TRUE`, temporarily raises the package's log level so
  `INFO`-level progress messages print for the duration of this call
  (default `FALSE` - quiet). See
  [`set_verbose_logging()`](https://jjmaynard.github.io/soilSIM/reference/set_verbose_logging.md).

## Value

Within-depth correlation validation
