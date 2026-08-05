# Fit Individual GP Model

Enhanced version with better error handling and diagnostics

## Usage

``` r
fit_individual_gp_model(
  data,
  property,
  optimize_hyperparameters = TRUE,
  verbose = getOption("ssurgo.verbose", FALSE)
)
```

## Arguments

- data:

  Input data for the group

- property:

  Property name to model

- optimize_hyperparameters:

  Whether to optimize hyperparameters

- verbose:

  Logical; if `TRUE`, temporarily raises the package's log level so
  `INFO`-level progress messages print for the duration of this call
  (default `FALSE` - quiet). See
  [`set_verbose_logging()`](https://jjmaynard.github.io/soilSIM/reference/set_verbose_logging.md).

## Value

List containing GP model and diagnostics
