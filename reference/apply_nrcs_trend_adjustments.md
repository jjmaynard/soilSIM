# Apply NRCS Trend Adjustments

Enhanced version with Module 8 integration and proper Module 5 function
calls.

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

  NRCS GP models from gp_modeling module

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
  `config$monte_carlo$vertical_correlation_method` (default
  `"joint_copula"` as of `VERTICAL_CORRELATION_IMPROVEMENT_PLAN.md`
  Phase 13; set to `"gp_quantile_retrofit"` to opt back into the
  original algorithm) reaches the NRCS/regional GP path the same way it
  already reaches the local-GP path (Phase 6/11). `NULL` (default)
  resolves to `"joint_copula"`, matching
  [`get_monte_carlo_defaults()`](https://jjmaynard.github.io/soilSIM/reference/get_monte_carlo_defaults.md)'s
  own default.

## Value

Adjusted simulation data
