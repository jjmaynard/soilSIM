# Apply NRCS Trend Adjustments

Maps each requested property to its NRCS regional GP model, predicts
that model's depth trend for the component's depths, and applies it via
apply_gp_depth_trends().

## Usage

``` r
apply_nrcs_trend_adjustments(
  cokey_data,
  gp_models,
  model_group,
  properties,
  preserve_correlations = TRUE,
  verbose = getOption("ssurgo.verbose", FALSE),
  config = NULL
)
```

## Arguments

- cokey_data:

  Simulation data for a single cokey

- gp_models:

  Fitted NRCS GP depth models

- model_group:

  GP model group for this cokey

- properties:

  Properties to adjust

- preserve_correlations:

  Whether to preserve correlations

- verbose:

  Logical; if `TRUE`, temporarily raises the package's log level so
  `INFO`-level progress messages print for the duration of this call
  (default `FALSE` - quiet). See
  [`set_verbose_logging()`](https://jjmaynard.github.io/soilSIM/reference/set_verbose_logging.md).

- config:

  Optional Monte Carlo config, passed through to
  [`apply_gp_depth_trends()`](https://jjmaynard.github.io/soilSIM/reference/apply_gp_depth_trends.md) -
  `config$monte_carlo$vertical_correlation_method` (`"joint_copula"`
  default, or `"gp_quantile_retrofit"`) reaches the NRCS/regional GP
  path as well as the local-GP path. `NULL` (default) resolves to
  `"joint_copula"`, matching
  [`get_monte_carlo_defaults()`](https://jjmaynard.github.io/soilSIM/reference/get_monte_carlo_defaults.md).

## Value

Adjusted simulation data
