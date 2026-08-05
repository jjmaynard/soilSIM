# Safe Correlation

Numerically stable correlation calculation with handling for edge cases.

## Usage

``` r
safe_correlation(
  x,
  y = NULL,
  method = "pearson",
  handle_constant = "warn",
  verbose = getOption("ssurgo.verbose", FALSE)
)
```

## Arguments

- x:

  First variable

- y:

  Second variable (NULL for correlation matrix)

- method:

  Correlation method ("pearson", "spearman", "kendall")

- handle_constant:

  How to handle constant variables ("remove", "zero", "warn")

- verbose:

  Logical; if `TRUE`, temporarily raises the package's log level so
  `INFO`-level progress messages print for the duration of this call
  (default `FALSE` - quiet). See
  [`set_verbose_logging()`](https://jjmaynard.github.io/soilSIM/reference/set_verbose_logging.md).

## Value

Correlation coefficient or matrix
