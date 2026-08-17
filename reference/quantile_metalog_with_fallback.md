# Metalog quantile with automatic fallback to `linear_cdf`

Metalog quantile with automatic fallback to `linear_cdf`

## Usage

``` r
quantile_metalog_with_fallback(fit, infeasible, full_probs, full_values, q)
```

## Arguments

- fit:

  Output of
  [`fit_metalog_linear()`](https://jjmaynard.github.io/soilSIM/reference/fit_metalog_linear.md).

- infeasible:

  Output of
  [`check_metalog_feasible()`](https://jjmaynard.github.io/soilSIM/reference/check_metalog_feasible.md)
  for this same `fit`.

- full_probs, full_values:

  The FULL percentile set (including p=0/p=1 if available) used for the
  `linear_cdf` fallback.

- q:

  Vector of probabilities.
