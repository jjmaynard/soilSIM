# Piecewise-linear inverse-CDF quantile function, exact at the given knots

Tails are clamped to the min/max value (`rule = 2`) rather than
extrapolated.

## Usage

``` r
quantile_linear_cdf(probs, values, q)
```

## Arguments

- probs, values:

  Matching, sorted percentile probabilities/values.

- q:

  Vector of probabilities to evaluate at.
