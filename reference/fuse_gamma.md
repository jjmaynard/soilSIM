# Fuse two independent Gamma (shape/rate) belief distributions

`Gamma(k1,r1) x Gamma(k2,r2) -> Gamma(k1+k2-1, r1+r2)`, derived from the
Gamma kernel `x^(k-1)exp(-r*x)`. Requires `k1+k2 > 1` to remain valid.

## Usage

``` r
fuse_gamma(prior_shape, prior_rate, lik_shape, lik_rate)
```

## Arguments

- prior_shape, prior_rate, lik_shape, lik_rate:

  Numeric scalars or vectors.

## Value

`list(shape=, rate=, feasible=)`.

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
[`fuse_normal_normal()`](https://jjmaynard.github.io/soilSIM/reference/fuse_normal_normal.md),
[`fuse_property()`](https://jjmaynard.github.io/soilSIM/reference/fuse_property.md),
[`update_prior()`](https://jjmaynard.github.io/soilSIM/reference/update_prior.md)
