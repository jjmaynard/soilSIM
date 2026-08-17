# Assess Simulation Convergence

Splits `simulation_data`'s realizations (`simulation_number`) into
sequential batches and checks whether each numeric property's running
mean has stabilized (relative change between the last two batches below
`criteria$tolerance`).

## Usage

``` r
assess_simulation_convergence(simulation_data, criteria)
```

## Arguments

- simulation_data:

  Monte Carlo simulation data (must have `simulation_number` and numeric
  property columns).

- criteria:

  List, optionally with `tolerance` (default 0.05) and `n_batches`
  (default 5).

## Value

List with `converged`, `convergence_metric` (max relative change across
assessed properties), `assessment_method`.
