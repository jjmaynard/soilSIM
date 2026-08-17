# Calculate Simulation Quality Metrics

Computes real per-property finite-value and outlier rates from the
`[horizon, property, realization]` simulation array, and an overall
quality score combining them.

## Usage

``` r
calculate_simulation_quality_metrics(simulation_results, properties, config)
```

## Arguments

- simulation_results:

  `[horizon, property, realization]` array.

- properties:

  Character vector of property names (column order of
  `simulation_results`'s second dimension).

- config:

  Simulation configuration (uses
  `config$monte_carlo$outlier_threshold`).

## Value

List with `per_property` (named list of `finite_rate`/`outlier_rate`),
`overall_finite_rate`, `overall_outlier_rate`, `overall_quality`.
