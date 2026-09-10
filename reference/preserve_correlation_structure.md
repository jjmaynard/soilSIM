# Preserve Correlation Structure During GP Adjustment

Imposes vertical correlation on the property matrices by re-ordering
each property's per-depth values to a shared surface-referenced rank, so
the depth trend is applied without disturbing within-depth
cross-property correlation.

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
