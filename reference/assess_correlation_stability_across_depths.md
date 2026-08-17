# Assess Correlation Stability Across Depths

Computes a correlation matrix for `numeric_properties` within each depth
bin of `simulation_data`, and measures how stable pairwise correlations
are across bins (1 - mean standard deviation of each pairwise
correlation across bins).

## Usage

``` r
assess_correlation_stability_across_depths(
  simulation_data,
  numeric_properties,
  criteria
)
```

## Arguments

- simulation_data:

  Simulation data with an `hzdept_r` column.

- numeric_properties:

  Character vector of property columns to correlate.

- criteria:

  List, optionally with `depth_bins` (default 5 quantile bins).

## Value

List with `stability_score`, `unstable_depths`, `stability_assessment`.
