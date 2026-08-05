# Apply Local GP Adjustments With Correlation Preservation

Delegates to the real
[`apply_local_gp_adjustments()`](https://jjmaynard.github.io/soilSIM/reference/apply_local_gp_adjustments.md)
(`multivariate-adjustment.R`, migrated from mod07), inferring the
`properties` argument that function requires (as the numeric,
non-structural columns of `simulation_data`), since this function's own
signature does not accept one explicitly.

## Usage

``` r
apply_local_gp_adjustments_with_correlations(
  simulation_data,
  preserve_correlations = TRUE
)
```

## Arguments

- simulation_data:

  Simulated soil property data for one cokey.

- preserve_correlations:

  Whether to preserve within-depth correlations.

## Value

Adjusted simulation data.
