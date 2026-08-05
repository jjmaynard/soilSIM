# Validate GP Models

Enhanced validation with comprehensive diagnostics

## Usage

``` r
validate_gp_models(
  gp_models,
  nrcs_data,
  validation_depths = seq(0, 200, by = 10),
  verbose = getOption("ssurgo.verbose", FALSE)
)
```

## Arguments

- gp_models:

  List of GP models from build_stratified_gp_models()

- nrcs_data:

  Original NRCS training data

- validation_depths:

  Depths for prediction testing

- verbose:

  Logical; if `TRUE`, temporarily raises the package's log level so
  `INFO`-level progress messages print for the duration of this call
  (default `FALSE` - quiet). See
  [`set_verbose_logging()`](https://jjmaynard.github.io/soilSIM/reference/set_verbose_logging.md).

## Value

Validation results
