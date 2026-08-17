# Piecewise-linear inverse-CDF quantile function, exact at the given knots

Tails are clamped to the min/max value (`rule = 2`) rather than
extrapolated, matching `code_ref/brdf/distribution_fitting.R`'s
[`sim_linear_cdf()`](https://jjmaynard.github.io/soilSIM/reference/sim_linear_cdf.md).

## Usage

``` r
quantile_linear_cdf(probs, values, q)
```

## Arguments

- probs, values:

  Matching, sorted percentile probabilities/values.

- q:

  Vector of probabilities to evaluate at.
