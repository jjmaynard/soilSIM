# Generate Simulation Diagnostics

Computes real per-property simulated vs. original mean/SD comparisons.

## Usage

``` r
generate_simulation_diagnostics(
  original_data,
  simulation_results,
  properties,
  correlation_config,
  config
)
```

## Arguments

- original_data:

  Original SSURGO-derived data.

- simulation_results:

  Constrained `[horizon, property, realization]` array.

- properties:

  Character vector of property names.

- correlation_config:

  Correlation configuration used for simulation.

- config:

  Unused (kept for interface compatibility).

## Value

List with `per_property` (named list of `sim_mean`/`sim_sd`/
`orig_mean`/`orig_sd`), `correlation_method`, `n_properties`,
`timestamp`.
