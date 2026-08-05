# Assess Component Composition Quality

Compares each component's mean simulated value against its original
`comppct_r`, as a relative-error-based quality score.

## Usage

``` r
assess_component_quality(original_data, simulated_data, config)
```

## Arguments

- original_data:

  Original component data with a `comppct_r` column.

- simulated_data:

  `[n_components, n_realizations]` matrix of simulated/normalized
  component values.

- config:

  Unused (kept for interface compatibility).

## Value

List with `overall_quality`, `mean_relative_error`.
