# Validate GP Model Workflow

Enhanced GP validation using Module 5 functions and Module 0 utilities.

## Usage

``` r
validate_gp_model_workflow(
  gp_models,
  training_data,
  config = NULL,
  verbose = getOption("ssurgo.verbose", FALSE)
)
```

## Arguments

- gp_models:

  GP models from gp_modeling module

- training_data:

  Original training data

- config:

  GP validation configuration

- verbose:

  Logical; if `TRUE`, temporarily raises the package's log level so
  `INFO`-level progress messages print for the duration of this call
  (default `FALSE` - quiet). See
  [`set_verbose_logging()`](https://jjmaynard.github.io/soilSIM/reference/set_verbose_logging.md).

## Value

GP model validation results
