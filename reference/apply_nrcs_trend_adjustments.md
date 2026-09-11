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
  [`default_monte_carlo_config()`](https://jjmaynard.github.io/soilSIM/reference/default_monte_carlo_config.md).

## Value

Adjusted simulation data

## See also

Other gp-modeling:
[`adjust_simulation_depthwise()`](https://jjmaynard.github.io/soilSIM/reference/adjust_simulation_depthwise.md),
[`apply_gp_depth_trends()`](https://jjmaynard.github.io/soilSIM/reference/apply_gp_depth_trends.md),
[`apply_local_gp_adjustments()`](https://jjmaynard.github.io/soilSIM/reference/apply_local_gp_adjustments.md),
[`convert_to_property_matrices()`](https://jjmaynard.github.io/soilSIM/reference/convert_to_property_matrices.md),
[`extract_depth_length_scale()`](https://jjmaynard.github.io/soilSIM/reference/extract_depth_length_scale.md),
[`fit_depth_gp_models()`](https://jjmaynard.github.io/soilSIM/reference/fit_depth_gp_models.md),
[`fit_individual_gp_model()`](https://jjmaynard.github.io/soilSIM/reference/fit_individual_gp_model.md),
[`fit_local_gp_model_single()`](https://jjmaynard.github.io/soilSIM/reference/fit_local_gp_model_single.md),
[`predict.soilSIM_gp_models()`](https://jjmaynard.github.io/soilSIM/reference/predict.soilSIM_gp_models.md),
[`predict_gp_depth_trends()`](https://jjmaynard.github.io/soilSIM/reference/predict_gp_depth_trends.md),
[`prepare_nrcs_training_data()`](https://jjmaynard.github.io/soilSIM/reference/prepare_nrcs_training_data.md),
[`select_optimal_grouping()`](https://jjmaynard.github.io/soilSIM/reference/select_optimal_grouping.md),
[`validate_gp_models()`](https://jjmaynard.github.io/soilSIM/reference/validate_gp_models.md)
