# Validate a Single GP Model's Predictions

Predicts a depth trend from `group_model` and measures its smoothness as
one minus the normalized mean absolute second difference of the
predicted curve (a roughness measure).

## Usage

``` r
validate_single_gp_predictions(group_model, criteria)
```

## Arguments

- group_model:

  Fitted GP model.

- criteria:

  List, optionally with `test_depths` (default `seq(0, 150, by = 10)`)
  and `smoothness_criteria` (default 0.1).

## Value

List with `prediction_smoothness`, `uncertainty_reasonable`,
`prediction_quality`.
