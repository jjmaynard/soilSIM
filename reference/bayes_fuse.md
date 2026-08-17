# Fuse a prior and likelihood belief distribution of the same family

Fuse a prior and likelihood belief distribution of the same family

## Usage

``` r
bayes_fuse(prior_params, lik_params, family = c("normal", "beta", "gamma"))
```

## Arguments

- prior_params, lik_params:

  Named lists of that family's parameters: `list(mu=,sigma=)` for
  `"normal"`, `list(alpha=,beta=)` for `"beta"`, `list(shape=,rate=)`
  for `"gamma"`.

- family:

  One of `"normal"`, `"beta"`, `"gamma"`. Both sides must already be fit
  in this same family - there is no closed form for fusing mismatched
  families; use
  [`bayesian_update()`](https://jjmaynard.github.io/soilSIM/reference/bayesian_update.md)
  for that.

## Value

The matching `fuse_*()` function's return value.
