# Validate GP Model Performance

Validate GP Model Performance

## Usage

``` r
diagnose_gp_performance(
  gp_models,
  training_data,
  criteria = NULL,
  verbose = getOption("ssurgo.verbose", FALSE)
)
```

## Arguments

- gp_models:

  GP models

- training_data:

  Training data

- criteria:

  Performance criteria

- verbose:

  Logical; if `TRUE`, temporarily raises the package's log level so
  `INFO`-level progress messages print for the duration of this call
  (default `FALSE` - quiet). See
  [`set_verbose_logging()`](https://jjmaynard.github.io/soilSIM/reference/set_verbose_logging.md).

## Value

GP model performance assessment
