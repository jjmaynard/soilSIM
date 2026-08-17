# Compare Distribution Moments

Compares `simulated_values`' mean/variance against a reference sample
drawn from the fitted distribution.

## Usage

``` r
compare_distribution_moments(simulated_values, param_info, criteria)
```

## Arguments

- simulated_values:

  Simulated values.

- param_info:

  Fitted distribution parameters, with `family` and `fit`.

- criteria:

  List, optionally with `moment_tolerance` (default 0.1).

## Value

List with `mean_difference`, `variance_difference` (relative),
`moment_match_quality`.
