# Inverse of `normal_to_lognormal_params()`

Converts log-space Normal(mu, sigma) parameters back to the raw-space
Lognormal's mean/sd.

## Usage

``` r
lognormal_to_normal_params(mu_log, sigma_log)
```

## Arguments

- mu_log, sigma_log:

  Log-space mean/sd.

## Value

`list(mu = raw-space mean, sigma = raw-space sd)`.
