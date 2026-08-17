# Predict GP Depth Trends

Enhanced version with better error handling

## Usage

``` r
predict_gp_depth_trends(
  gp_model_info,
  new_depths,
  verbose = getOption("ssurgo.verbose", FALSE)
)
```

## Arguments

- gp_model_info:

  GP model information from fit_individual_gp_model()

- new_depths:

  Vector of depths for prediction

- verbose:

  Logical; if `TRUE`, temporarily raises the package's log level so
  `INFO`-level progress messages print for the duration of this call
  (default `FALSE` - quiet). See
  [`set_verbose_logging()`](https://jjmaynard.github.io/soilSIM/reference/set_verbose_logging.md).

## Value

Vector of predictions
