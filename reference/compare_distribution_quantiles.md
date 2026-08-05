# Compare Distribution Quantiles

Compares `simulated_values`' 10th/50th/90th percentiles against the
fitted distribution's theoretical quantiles at those same probabilities.

## Usage

``` r
compare_distribution_quantiles(simulated_values, param_info, criteria)
```

## Arguments

- simulated_values:

  Simulated values.

- param_info:

  Fitted distribution parameters, with `family` and `fit`.

- criteria:

  List, optionally with `quantile_tolerance` (default 0.05).

## Value

List with `quantile_differences` (relative, length 3),
`quantile_match_quality`.
