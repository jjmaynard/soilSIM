# Run Monte Carlo simulation across multiple worker processes

simulate_correlated_properties() generates all n_realizations for every
horizon in one call, with no state carried across horizons or
realizations, so it can safely be called once per worker with a smaller
n_realizations chunk and the resulting arrays concatenated along the
realization dimension. Follows the same Windows-PSOCK-vs-mclapply
pattern as mod07_multivariate_adjustment.R's process_cokeys_parallel(),
with a sequential fallback on any error.

## Usage

``` r
run_parallel_simulation(
  simulation_params,
  distribution_setup,
  correlation_config,
  n_realizations,
  properties,
  config
)
```

## Arguments

- simulation_params:

  List of per-horizon parameter lists.

- distribution_setup:

  Unused; kept for interface compatibility with
  `run_sequential_simulation()`.

- correlation_config:

  Result of
  [`configure_correlation_structure()`](https://jjmaynard.github.io/soilSIM/reference/configure_correlation_structure.md).

- n_realizations:

  Number of realizations to generate.

- properties:

  Character vector of property names.

- config:

  Simulation configuration.

## Value

Same shape as simulate_correlated_properties(): an array dimensioned
horizon by property by realization.
