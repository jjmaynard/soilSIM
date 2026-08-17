# Validate Correlation Preservation

Enhanced correlation preservation validation using Module 6 utilities.

## Usage

``` r
validate_correlation_preservation_diagnostics(
  original_correlations,
  simulation_data,
  criteria = NULL,
  verbose = getOption("ssurgo.verbose", FALSE)
)
```

## Arguments

- original_correlations:

  Original correlation matrices

- simulation_data:

  Final simulation data

- criteria:

  Preservation criteria

- verbose:

  Logical; if `TRUE`, temporarily raises the package's log level so
  `INFO`-level progress messages print for the duration of this call
  (default `FALSE` - quiet). See
  [`set_verbose_logging()`](https://jjmaynard.github.io/soilSIM/reference/set_verbose_logging.md).

## Value

Correlation preservation assessment
