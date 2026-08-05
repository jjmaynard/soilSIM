# Validate Simulation Depth Trends

Compares depth-binned means of `simulation_data` against `original_data`
for shared numeric properties, via the correlation between simulated and
original per-depth means.

## Usage

``` r
validate_simulation_depth_trends(simulation_data, original_data, criteria)
```

## Arguments

- simulation_data:

  Simulation data.

- original_data:

  Original data for comparison.

- criteria:

  Unused.

## Value

List with `depth_trend_realism` (mean rescaled correlation across
properties), `trend_violations` (count of properties with correlation
below 0.5), `overall_trend_quality`.
