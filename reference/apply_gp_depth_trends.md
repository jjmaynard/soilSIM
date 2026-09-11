# Apply GP Depth Trends with Correlation Preservation

Applies GP-predicted depth trends to one component's simulation data,
optionally preserving within-depth cross-property correlations.

## Usage

``` r
apply_gp_depth_trends(
  cokey_data,
  gp_predictions,
  properties,
  preserve_correlations = TRUE,
  primary_property = NULL,
  verbose = getOption("ssurgo.verbose", FALSE),
  config = NULL,
  gp_models = NULL
)
```

## Arguments

- cokey_data:

  Simulation data for a single cokey

- gp_predictions:

  Named list of GP predictions by property

- properties:

  Properties to adjust

- preserve_correlations:

  Whether to preserve correlations

- primary_property:

  Reference property for correlation preservation

- verbose:

  Logical; if `TRUE`, temporarily raises the package's log level so
  `INFO`-level progress messages print for the duration of this call
  (default `FALSE` - quiet). See
  [`set_verbose_logging()`](https://jjmaynard.github.io/soilSIM/reference/set_verbose_logging.md).

- config:

  Optional Monte Carlo config (as from
  [`default_monte_carlo_config()`](https://jjmaynard.github.io/soilSIM/reference/default_monte_carlo_config.md))
  whose `monte_carlo$vertical_correlation_method` selects between
  `"joint_copula"` (the default; dispatches to
  [`preserve_correlation_structure_joint()`](https://jjmaynard.github.io/soilSIM/reference/preserve_correlation_structure_joint.md),
  drawing depth correlation and property correlation simultaneously) and
  `"gp_quantile_retrofit"` (an explicit opt-out that dispatches to
  [`preserve_correlation_structure()`](https://jjmaynard.github.io/soilSIM/reference/preserve_correlation_structure.md)).
  `NULL` (default) resolves to `"joint_copula"`, matching
  [`default_monte_carlo_config()`](https://jjmaynard.github.io/soilSIM/reference/default_monte_carlo_config.md)'s
  own default - set
  `config$monte_carlo$vertical_correlation_method = "gp_quantile_retrofit"`
  explicitly to select the quantile-retrofit method. Under
  `"joint_copula"`, `config$monte_carlo$vertical_correlation_gating`
  (default `FALSE`) separately controls whether `bound_sd`-based
  discontinuity gating is applied. Its numeric defaults are conservative
  and kept independent of the core method choice.

- gp_models:

  Optional named list of fitted GP models (as
  [`fit_local_gp_model_single()`](https://jjmaynard.github.io/soilSIM/reference/fit_local_gp_model_single.md)
  returns), keyed by property - passed through to
  [`preserve_correlation_structure_joint()`](https://jjmaynard.github.io/soilSIM/reference/preserve_correlation_structure_joint.md)
  when `vertical_correlation_method = "joint_copula"`, so its depth
  kernel can reuse each property's already-fitted,
  already-cross-validated length-scale (see
  [`extract_depth_length_scale()`](https://jjmaynard.github.io/soilSIM/reference/extract_depth_length_scale.md))
  instead of falling back to a full-depth-range default. Ignored
  entirely under `"gp_quantile_retrofit"`.

## Value

Adjusted simulation data

## See also

Other gp-modeling:
[`adjust_simulation_depthwise()`](https://jjmaynard.github.io/soilSIM/reference/adjust_simulation_depthwise.md),
[`apply_local_gp_adjustments()`](https://jjmaynard.github.io/soilSIM/reference/apply_local_gp_adjustments.md),
[`apply_nrcs_trend_adjustments()`](https://jjmaynard.github.io/soilSIM/reference/apply_nrcs_trend_adjustments.md),
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
