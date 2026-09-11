# Predict a fitted GP depth-trend model at new depths

A convenience wrapper over
[`predict_gp_depth_trends()`](https://jjmaynard.github.io/soilSIM/reference/predict_gp_depth_trends.md)
(unchanged - this does not replace it, since that function's real
callers throughout the package operate on one leaf model at a time, not
the whole `soilSIM_gp_models` collection) that looks up the requested
property/group's fitted model inside a
[`fit_depth_gp_models()`](https://jjmaynard.github.io/soilSIM/reference/fit_depth_gp_models.md)
result.

## Usage

``` r
# S3 method for class 'soilSIM_gp_models'
predict(object, new_depths, property, group = NULL, ...)
```

## Arguments

- object:

  A `soilSIM_gp_models` object from
  [`fit_depth_gp_models()`](https://jjmaynard.github.io/soilSIM/reference/fit_depth_gp_models.md).

- new_depths:

  Numeric vector of depths (cm) to predict at.

- property:

  Which property's model to use (a name in `object`).

- group:

  Which stratified group's model to use (default: the first group fitted
  for `property`, with a message when more than one exists - pass
  explicitly to silence it).

- ...:

  Passed to
  [`predict_gp_depth_trends()`](https://jjmaynard.github.io/soilSIM/reference/predict_gp_depth_trends.md).

## Value

Numeric vector of predicted means, one per `new_depths` (see
[`predict_gp_depth_trends()`](https://jjmaynard.github.io/soilSIM/reference/predict_gp_depth_trends.md)'s
own `@return` - it is the fitted mean only, not a prediction interval).

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
[`predict_gp_depth_trends()`](https://jjmaynard.github.io/soilSIM/reference/predict_gp_depth_trends.md),
[`prepare_nrcs_training_data()`](https://jjmaynard.github.io/soilSIM/reference/prepare_nrcs_training_data.md),
[`select_optimal_grouping()`](https://jjmaynard.github.io/soilSIM/reference/select_optimal_grouping.md),
[`validate_gp_models()`](https://jjmaynard.github.io/soilSIM/reference/validate_gp_models.md)
