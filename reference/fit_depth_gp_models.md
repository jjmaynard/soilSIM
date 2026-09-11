# Build Stratified GP Models

Build Stratified GP Models

## Usage

``` r
fit_depth_gp_models(
  processed_nrcs_data,
  properties = c("clay_pct", "sand_pct", "pH", "organic_matter"),
  min_profiles_per_group = 3,
  min_observations_per_group = 15,
  optimize_hyperparameters = TRUE,
  verbose = getOption("ssurgo.verbose", FALSE)
)
```

## Arguments

- processed_nrcs_data:

  Processed NRCS data from prepare_nrcs_training_data()

- properties:

  Character vector of properties to model

- min_profiles_per_group:

  Minimum profiles per group (default = 3)

- min_observations_per_group:

  Minimum observations per group (default = 15)

- optimize_hyperparameters:

  Whether to optimize GP hyperparameters (default = TRUE)

- verbose:

  Logical; if `TRUE`, temporarily raises the package's log level so
  `INFO`-level progress messages print for the duration of this call
  (default `FALSE` - quiet). See
  [`set_verbose_logging()`](https://jjmaynard.github.io/soilSIM/reference/set_verbose_logging.md).

## Value

List of GP models organized by property and group

## See also

Other gp-modeling:
[`adjust_simulation_depthwise()`](https://jjmaynard.github.io/soilSIM/reference/adjust_simulation_depthwise.md),
[`apply_gp_depth_trends()`](https://jjmaynard.github.io/soilSIM/reference/apply_gp_depth_trends.md),
[`apply_local_gp_adjustments()`](https://jjmaynard.github.io/soilSIM/reference/apply_local_gp_adjustments.md),
[`apply_nrcs_trend_adjustments()`](https://jjmaynard.github.io/soilSIM/reference/apply_nrcs_trend_adjustments.md),
[`convert_to_property_matrices()`](https://jjmaynard.github.io/soilSIM/reference/convert_to_property_matrices.md),
[`extract_depth_length_scale()`](https://jjmaynard.github.io/soilSIM/reference/extract_depth_length_scale.md),
[`fit_individual_gp_model()`](https://jjmaynard.github.io/soilSIM/reference/fit_individual_gp_model.md),
[`fit_local_gp_model_single()`](https://jjmaynard.github.io/soilSIM/reference/fit_local_gp_model_single.md),
[`predict.soilSIM_gp_models()`](https://jjmaynard.github.io/soilSIM/reference/predict.soilSIM_gp_models.md),
[`predict_gp_depth_trends()`](https://jjmaynard.github.io/soilSIM/reference/predict_gp_depth_trends.md),
[`prepare_nrcs_training_data()`](https://jjmaynard.github.io/soilSIM/reference/prepare_nrcs_training_data.md),
[`select_optimal_grouping()`](https://jjmaynard.github.io/soilSIM/reference/select_optimal_grouping.md),
[`validate_gp_models()`](https://jjmaynard.github.io/soilSIM/reference/validate_gp_models.md)
