# Assess Property Coverage

Compares the range/quantile coverage of `sim_values` against
`orig_values`: what fraction of `orig_values`' range falls within
`sim_values`' range, and whether `sim_values` extrapolates meaningfully
beyond it.

## Usage

``` r
assess_property_coverage(sim_values, orig_values, criteria)
```

## Arguments

- sim_values:

  Simulated values.

- orig_values:

  Original/reference values.

- criteria:

  List, optionally with `max_extrapolation_factor` (default 1.2).

## Value

List with `coverage_percentage`, `extrapolation_detected`,
`coverage_quality`.
