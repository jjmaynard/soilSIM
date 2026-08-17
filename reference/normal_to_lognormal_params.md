# Moment-match a raw-space Lognormal's mean/sd onto its underlying Normal's mu/sigma

Needed because several properties (e.g. SOC, CEC) are strictly positive
and right-skewed - fusing their raw mean/sd as if Normal lets the
posterior put real probability mass below zero, which is physically
impossible.
[`bayes_update_normal_normal()`](https://jjmaynard.github.io/soilSIM/reference/bayes_update_normal_normal.md)
is reused in log-space instead.

## Usage

``` r
normal_to_lognormal_params(mu, sigma)
```

## Arguments

- mu, sigma:

  Raw-space mean/sd.

## Value

`list(mu = log-space mu, sigma = log-space sigma)`.
