# Prepare NRCS Training Data for GP Model Building

Prepare NRCS Training Data for GP Model Building

## Usage

``` r
prepare_nrcs_training_data(
  nrcs_combined_data,
  grouping_strategy = "auto",
  min_profiles_per_group = 3,
  min_observations_per_group = 15,
  target_min_groups = 3,
  max_depth = DEFAULT_MAX_DEPTH_CM,
  validation_config = NULL,
  verbose = getOption("ssurgo.verbose", FALSE)
)
```

## Arguments

- nrcs_combined_data:

  Combined NRCS data with horizon and component information

- grouping_strategy:

  How to group soils ("auto" for intelligent selection)

- min_profiles_per_group:

  Minimum profiles required per group (default = 3)

- min_observations_per_group:

  Minimum total observations per group (default = 15)

- target_min_groups:

  Minimum number of adequate groups desired (default = 3)

- max_depth:

  Maximum depth to include in training (default = 250 cm)

- validation_config:

  Validation configuration from the shared utilities

- verbose:

  Logical; if `TRUE`, temporarily raises the package's log level so
  `INFO`-level progress messages print for the duration of this call
  (default `FALSE` - quiet). See
  [`set_verbose_logging()`](https://jjmaynard.github.io/soilSIM/reference/set_verbose_logging.md).

## Value

Processed data frame ready for GP model building

## See also

Other gp-modeling:
[`adjust_simulation_depthwise()`](https://jjmaynard.github.io/soilSIM/reference/adjust_simulation_depthwise.md),
[`apply_gp_depth_trends()`](https://jjmaynard.github.io/soilSIM/reference/apply_gp_depth_trends.md),
[`apply_local_gp_adjustments()`](https://jjmaynard.github.io/soilSIM/reference/apply_local_gp_adjustments.md),
[`apply_nrcs_trend_adjustments()`](https://jjmaynard.github.io/soilSIM/reference/apply_nrcs_trend_adjustments.md),
[`convert_to_property_matrices()`](https://jjmaynard.github.io/soilSIM/reference/convert_to_property_matrices.md),
[`extract_depth_length_scale()`](https://jjmaynard.github.io/soilSIM/reference/extract_depth_length_scale.md),
[`fit_depth_gp_models()`](https://jjmaynard.github.io/soilSIM/reference/fit_depth_gp_models.md),
[`fit_individual_gp_model()`](https://jjmaynard.github.io/soilSIM/reference/fit_individual_gp_model.md),
[`fit_local_gp_model_single()`](https://jjmaynard.github.io/soilSIM/reference/fit_local_gp_model_single.md),
[`predict.soilSIM_gp_models()`](https://jjmaynard.github.io/soilSIM/reference/predict.soilSIM_gp_models.md),
[`predict_gp_depth_trends()`](https://jjmaynard.github.io/soilSIM/reference/predict_gp_depth_trends.md),
[`select_optimal_grouping()`](https://jjmaynard.github.io/soilSIM/reference/select_optimal_grouping.md),
[`validate_gp_models()`](https://jjmaynard.github.io/soilSIM/reference/validate_gp_models.md)
