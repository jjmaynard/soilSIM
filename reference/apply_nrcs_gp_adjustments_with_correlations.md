# Apply NRCS GP Adjustments With Correlation Preservation

Delegates to the real
[`apply_nrcs_trend_adjustments()`](https://jjmaynard.github.io/soilSIM/reference/apply_nrcs_trend_adjustments.md)
(`multivariate-adjustment.R`), which already implements NRCS GP-model
depth-trend adjustment with correlation preservation. This wrapper's
only job is to infer the `properties` argument that function requires
(as the numeric, non-structural columns of `simulation_data`), since
this function's own signature does not accept one explicitly.

## Usage

``` r
apply_nrcs_gp_adjustments_with_correlations(
  simulation_data,
  nrcs_gp_models,
  model_group,
  preserve_correlations = TRUE
)
```

## Arguments

- simulation_data:

  Simulated soil property data for one cokey.

- nrcs_gp_models:

  Fitted NRCS GP models (see
  [`build_stratified_gp_models()`](https://jjmaynard.github.io/soilSIM/reference/build_stratified_gp_models.md)).

- model_group:

  GP model group to apply.

- preserve_correlations:

  Whether to preserve within-depth correlations.

## Value

Adjusted simulation data.
