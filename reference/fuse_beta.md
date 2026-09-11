# Fuse two independent Beta belief distributions via density multiplication

`Beta(a1,b1) x Beta(a2,b2) -> Beta(a1+a2-1, b1+b2-1)`, derived from the
Beta kernel `x^(a-1)(1-x)^(b-1)`. Requires `a1+a2 > 1` and `b1+b2 > 1`
to remain a valid distribution; `feasible` flags where that fails.

## Usage

``` r
fuse_beta(prior_alpha, prior_beta, lik_alpha, lik_beta)
```

## Arguments

- prior_alpha, prior_beta, lik_alpha, lik_beta:

  Numeric scalars or vectors.

## Value

`list(alpha=, beta=, feasible=)`.

## See also

Other bayesian-fusion:
[`convert_beta_to_moments()`](https://jjmaynard.github.io/soilSIM/reference/convert_beta_to_moments.md),
[`convert_gamma_to_moments()`](https://jjmaynard.github.io/soilSIM/reference/convert_gamma_to_moments.md),
[`convert_lognormal_to_normal()`](https://jjmaynard.github.io/soilSIM/reference/convert_lognormal_to_normal.md),
[`convert_moments_to_beta()`](https://jjmaynard.github.io/soilSIM/reference/convert_moments_to_beta.md),
[`convert_moments_to_gamma()`](https://jjmaynard.github.io/soilSIM/reference/convert_moments_to_gamma.md),
[`convert_normal_to_lognormal()`](https://jjmaynard.github.io/soilSIM/reference/convert_normal_to_lognormal.md),
[`fuse_bivariate_normal()`](https://jjmaynard.github.io/soilSIM/reference/fuse_bivariate_normal.md),
[`fuse_distribution()`](https://jjmaynard.github.io/soilSIM/reference/fuse_distribution.md),
[`fuse_gamma()`](https://jjmaynard.github.io/soilSIM/reference/fuse_gamma.md),
[`fuse_normal_normal()`](https://jjmaynard.github.io/soilSIM/reference/fuse_normal_normal.md),
[`fuse_property()`](https://jjmaynard.github.io/soilSIM/reference/fuse_property.md),
[`update_prior()`](https://jjmaynard.github.io/soilSIM/reference/update_prior.md)
