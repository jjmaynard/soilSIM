# Precision-weighted conjugate Normal-Normal posterior

Precision-weighted conjugate Normal-Normal posterior

## Usage

``` r
bayes_update_normal_normal(prior_mu, prior_sigma, lik_mu, lik_sigma)
```

## Arguments

- prior_mu, prior_sigma:

  Prior mean/sd (numeric scalar or vector).

- lik_mu, lik_sigma:

  Likelihood mean/sd (numeric scalar or vector, same length as the prior
  args).

## Value

`list(mu = posterior mean, sigma = posterior sd)`.

## Details

Invariants (see `tests/testthat/test-bayesian-updating.R`): posterior
`sigma <= min(prior_sigma, lik_sigma)`; posterior `mu` is bounded
between `prior_mu` and `lik_mu`; symmetric under swapping
`(prior_mu, prior_sigma)` \<-\> `(lik_mu, lik_sigma)`.
