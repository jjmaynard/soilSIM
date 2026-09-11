# Apply Local GP Adjustments

Fits per-component Gaussian-process depth models for the requested
properties and applies their predicted trends to the component's
simulation data.

## Usage

``` r
apply_local_gp_adjustments(
  cokey_data,
  properties,
  preserve_correlations = TRUE,
  min_depths = 3,
  config = NULL,
  gp_control = c(20, 10, 2),
  verbose = getOption("ssurgo.verbose", FALSE)
)
```

## Arguments

- cokey_data:

  Simulation data for a single cokey

- properties:

  Properties to adjust

- preserve_correlations:

  Whether to preserve correlations

- min_depths:

  Minimum depths required for GP fitting

- config:

  Optional configuration list; defaults to the package validation
  configuration

- gp_control:

  Passed through to
  [`fit_local_gp_models()`](https://jjmaynard.github.io/soilSIM/reference/fit_local_gp_models.md)/[`fit_local_gp_model_single()`](https://jjmaynard.github.io/soilSIM/reference/fit_local_gp_model_single.md)'s
  `gp_control` - see
  [`fit_local_gp_model_single()`](https://jjmaynard.github.io/soilSIM/reference/fit_local_gp_model_single.md)'s
  docs for why the default is much smaller than
  [`GPfit::GP_fit()`](https://rdrr.io/pkg/GPfit/man/GP_fit.html)'s own
  default.

- verbose:

  Logical; if `TRUE`, temporarily raises the package's log level so
  `INFO`-level progress messages print for the duration of this call
  (default `FALSE` - quiet). See
  [`set_verbose_logging()`](https://jjmaynard.github.io/soilSIM/reference/set_verbose_logging.md).

## Value

Adjusted simulation data

## See also

Other gp-modeling:
[`adjust_simulation_depthwise()`](https://jjmaynard.github.io/soilSIM/reference/adjust_simulation_depthwise.md),
[`apply_gp_depth_trends()`](https://jjmaynard.github.io/soilSIM/reference/apply_gp_depth_trends.md),
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
