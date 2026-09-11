# Fuse a prior and likelihood belief distribution of the same family

Fuse a prior and likelihood belief distribution of the same family

## Usage

``` r
fuse_distribution(
  prior_params,
  lik_params,
  family = c("normal", "beta", "gamma")
)
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
  [`update_prior()`](https://jjmaynard.github.io/soilSIM/reference/update_prior.md)
  for that.

## Value

The matching `fuse_*()` function's return value.

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
[`fuse_gamma()`](https://jjmaynard.github.io/soilSIM/reference/fuse_gamma.md),
[`fuse_normal_normal()`](https://jjmaynard.github.io/soilSIM/reference/fuse_normal_normal.md),
[`fuse_property()`](https://jjmaynard.github.io/soilSIM/reference/fuse_property.md),
[`update_prior()`](https://jjmaynard.github.io/soilSIM/reference/update_prior.md)
