# Precision-weighted conjugate Normal-Normal posterior

Precision-weighted conjugate Normal-Normal posterior

## Usage

``` r
fuse_normal_normal(prior_mu, prior_sigma, lik_mu, lik_sigma)
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

## See also

Other bayesian-fusion:
[`convert_beta_to_moments()`](https://jjmaynard.github.io/soilSIM/reference/convert_beta_to_moments.md),
[`convert_gamma_to_moments()`](https://jjmaynard.github.io/soilSIM/reference/convert_gamma_to_moments.md),
[`convert_lognormal_to_normal()`](https://jjmaynard.github.io/soilSIM/reference/convert_lognormal_to_normal.md),
[`convert_moments_to_beta()`](https://jjmaynard.github.io/soilSIM/reference/convert_moments_to_beta.md),
[`convert_moments_to_gamma()`](https://jjmaynard.github.io/soilSIM/reference/convert_moments_to_gamma.md),
[`convert_normal_to_lognormal()`](https://jjmaynard.github.io/soilSIM/reference/convert_normal_to_lognormal.md),
[`fuse_beta()`](https://jjmaynard.github.io/soilSIM/reference/fuse_beta.md),
[`fuse_bivariate_normal()`](https://jjmaynard.github.io/soilSIM/reference/fuse_bivariate_normal.md),
[`fuse_distribution()`](https://jjmaynard.github.io/soilSIM/reference/fuse_distribution.md),
[`fuse_gamma()`](https://jjmaynard.github.io/soilSIM/reference/fuse_gamma.md),
[`fuse_property()`](https://jjmaynard.github.io/soilSIM/reference/fuse_property.md),
[`update_prior()`](https://jjmaynard.github.io/soilSIM/reference/update_prior.md)
