# Adjust Multiple Soil Properties While Preserving Correlations

Core function from GP-depth-adjust.R that adjusts multiple simulated
soil properties simultaneously to follow GP-predicted trends while
preserving within-depth correlations.

## Usage

``` r
adjust_multivariate_depthwise_GP(
  simulated_list,
  gp_models,
  depths,
  primary_property = NULL,
  verbose = getOption("ssurgo.verbose", FALSE)
)
```

## Arguments

- simulated_list:

  A named list of matrices, where each matrix contains simulated values
  for one property (rows = depths, columns = simulations)

- gp_models:

  A named list of fitted GP models for each property

- depths:

  A numeric vector of depth values

- primary_property:

  Character string specifying which property to use as the "reference"
  for correlation preservation (default: first property in list)

- verbose:

  Logical; if `TRUE`, temporarily raises the package's log level so
  `INFO`-level progress messages print for the duration of this call
  (default `FALSE` - quiet). See
  [`set_verbose_logging()`](https://jjmaynard.github.io/soilSIM/reference/set_verbose_logging.md).

## Value

A named list of adjusted matrices with preserved correlations
