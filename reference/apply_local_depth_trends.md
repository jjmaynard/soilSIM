# Apply Local Depth Trends

Applies locally fitted GP depth-trend predictions to a component's
simulation data via apply_gp_depth_trends().

## Usage

``` r
apply_local_depth_trends(
  cokey_data,
  local_predictions,
  unique_depths,
  preserve_correlations = TRUE,
  verbose = getOption("ssurgo.verbose", FALSE),
  config = NULL,
  gp_models = NULL
)
```

## Arguments

- cokey_data:

  Simulation data

- local_predictions:

  Local GP predictions

- unique_depths:

  Depth vector

- preserve_correlations:

  Whether to preserve correlations

- verbose:

  Logical; if `TRUE`, temporarily raises the package's log level so
  `INFO`-level progress messages print for the duration of this call
  (default `FALSE` - quiet). See
  [`set_verbose_logging()`](https://jjmaynard.github.io/soilSIM/reference/set_verbose_logging.md).

- config:

  Optional config, passed straight through to
  [`apply_gp_depth_trends()`](https://jjmaynard.github.io/soilSIM/reference/apply_gp_depth_trends.md) -
  lets `config$monte_carlo$vertical_correlation_method`
  (`"joint_copula"` default, or `"gp_quantile_retrofit"`) reach this
  call site. `NULL` (default) resolves to `"joint_copula"`, matching
  [`get_monte_carlo_defaults()`](https://jjmaynard.github.io/soilSIM/reference/get_monte_carlo_defaults.md).

- gp_models:

  Optional named list of fitted local GP models (as
  [`apply_local_gp_adjustments()`](https://jjmaynard.github.io/soilSIM/reference/apply_local_gp_adjustments.md)
  already has in scope via
  [`fit_local_gp_models()`](https://jjmaynard.github.io/soilSIM/reference/fit_local_gp_models.md)),
  passed straight through to
  [`apply_gp_depth_trends()`](https://jjmaynard.github.io/soilSIM/reference/apply_gp_depth_trends.md)
  so the joint-copula depth kernel can reuse their fitted length-scales.
  Ignored under `"gp_quantile_retrofit"`.

## Value

Adjusted simulation data
