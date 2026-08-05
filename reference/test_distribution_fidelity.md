# Test Distribution Fidelity

Compares `simulated_values` against a reference sample drawn from the
fitted distribution (`param_info$family`/`param_info$fit`, via
[`quantile_from_fit()`](https://jjmaynard.github.io/soilSIM/reference/quantile_from_fit.md))
using a two-sample Kolmogorov-Smirnov test.

## Usage

``` r
test_distribution_fidelity(simulated_values, param_info, criteria)
```

## Arguments

- simulated_values:

  Simulated values.

- param_info:

  Fitted distribution parameters, with `family` and `fit`.

- criteria:

  List, optionally with `ks_test_alpha` (default 0.05).

## Value

List with `ks_test_p_value`, `distribution_match`, `test_method`.
