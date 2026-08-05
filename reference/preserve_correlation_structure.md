# Preserve Correlation Structure During GP Adjustment

Enhanced core correlation preservation algorithm with Module 8 error
handling.

## Usage

``` r
preserve_correlation_structure(
  property_matrices,
  gp_predictions,
  depths,
  primary_property,
  verbose = getOption("ssurgo.verbose", FALSE)
)
```

## Arguments

- property_matrices:

  Named list of property matrices

- gp_predictions:

  Named list of GP predictions

- depths:

  Depth vector

- primary_property:

  Reference property for correlation preservation

- verbose:

  Logical; if `TRUE`, temporarily raises the package's log level so
  `INFO`-level progress messages print for the duration of this call
  (default `FALSE` - quiet). See
  [`set_verbose_logging()`](https://jjmaynard.github.io/soilSIM/reference/set_verbose_logging.md).

## Value

List of adjusted property matrices
