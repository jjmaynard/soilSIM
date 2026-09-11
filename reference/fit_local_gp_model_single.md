# Fit a Single Local GP Model

Fit a Single Local GP Model

## Usage

``` r
fit_local_gp_model_single(agg_data, prop, gp_control = c(20, 10, 2))
```

## Arguments

- agg_data:

  One row per depth (`hzdept_r`, `mean_val`), as returned by
  `aggregate_property_by_depth()` - always a handful of points
  (typically \<10, one per unique depth in a single cokey), never one
  row per Monte Carlo replicate.

- prop:

  Property name, stored on the returned model for reference.

- gp_control:

  [`GPfit::GP_fit()`](https://rdrr.io/pkg/GPfit/man/GP_fit.html)'s
  `control` argument (population size / iteration counts for its
  internal hyperparameter search). Defaults to `c(20, 10, 2)`, far below
  `GP_fit()`'s own default `c(200*d, 80*d, 2*d)` (`d` = input
  dimensionality, always 1 here - depth is the only predictor).
  `GP_fit()`'s default search effort scales with `d`, not with the
  number of training points, and a 1-D, single-hyperparameter (`beta`)
  fit on this few points has a simple enough likelihood surface that the
  smaller search converges to the identical optimum every time -
  verified empirically across several representative depth/value series
  (identical fitted `beta` and identical predictions vs. the default,
  ~5x faster per call). Pass `c(200, 80, 2)` (or larger) to restore
  `GP_fit()`'s own default search effort if a future property/dataset
  needs a more thorough search.

## Value

A list with the fitted GP model, depth scaling info, training data, and
`prop`.

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
[`predict.soilSIM_gp_models()`](https://jjmaynard.github.io/soilSIM/reference/predict.soilSIM_gp_models.md),
[`predict_gp_depth_trends()`](https://jjmaynard.github.io/soilSIM/reference/predict_gp_depth_trends.md),
[`prepare_nrcs_training_data()`](https://jjmaynard.github.io/soilSIM/reference/prepare_nrcs_training_data.md),
[`select_optimal_grouping()`](https://jjmaynard.github.io/soilSIM/reference/select_optimal_grouping.md),
[`validate_gp_models()`](https://jjmaynard.github.io/soilSIM/reference/validate_gp_models.md)
