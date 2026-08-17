# Assess Cholesky Decomposition

Enhanced Cholesky validation error handling.

## Usage

``` r
assess_cholesky_decomposition(
  correlation_matrices,
  criteria = NULL,
  verbose = getOption("ssurgo.verbose", FALSE)
)
```

## Arguments

- correlation_matrices:

  Correlation matrices and decompositions

- criteria:

  Cholesky validation criteria

- verbose:

  Logical; if `TRUE`, temporarily raises the package's log level so
  `INFO`-level progress messages print for the duration of this call
  (default `FALSE` - quiet). See
  [`set_verbose_logging()`](https://jjmaynard.github.io/soilSIM/reference/set_verbose_logging.md).

## Value

Cholesky decomposition assessment
