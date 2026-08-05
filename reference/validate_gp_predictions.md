# Validate GP Predictions

Enhanced GP prediction validation using Module 5 functions.

## Usage

``` r
validate_gp_predictions(
  gp_models,
  criteria = NULL,
  verbose = getOption("ssurgo.verbose", FALSE)
)
```

## Arguments

- gp_models:

  GP models

- criteria:

  Prediction criteria

- verbose:

  Logical; if `TRUE`, temporarily raises the package's log level so
  `INFO`-level progress messages print for the duration of this call
  (default `FALSE` - quiet). See
  [`set_verbose_logging()`](https://jjmaynard.github.io/soilSIM/reference/set_verbose_logging.md).

## Value

GP prediction validation
