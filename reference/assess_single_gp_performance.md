# Assess a Single GP Model's Performance

Delegates to the already-real `calculate_model_diagnostics()`
(`gp-modeling.R`) for training RMSE, and derives an R-squared from that
RMSE against the training data's own variance.

## Usage

``` r
assess_single_gp_performance(group_model, criteria)
```

## Arguments

- group_model:

  Fitted GP model (as produced by
  [`fit_individual_gp_model()`](https://jjmaynard.github.io/soilSIM/reference/fit_individual_gp_model.md):
  `gp_model`, `training_data`, `property`).

- criteria:

  List, optionally with `max_training_rmse` (default `Inf`) and
  `min_r_squared` (default 0).

## Value

List with `training_rmse`, `r_squared`, `model_quality`,
`performance_assessment`.
