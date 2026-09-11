# Moment-match a raw-space Lognormal's mean/sd onto its underlying Normal's mu/sigma

Needed because several properties (e.g. SOC, CEC) are strictly positive
and right-skewed - fusing their raw mean/sd as if Normal lets the
posterior put real probability mass below zero, which is physically
impossible.
[`fuse_normal_normal()`](https://jjmaynard.github.io/soilSIM/reference/fuse_normal_normal.md)
is reused in log-space instead.

## Usage

``` r
convert_normal_to_lognormal(mu, sigma)
```

## Arguments

- mu, sigma:

  Raw-space mean/sd.

## Value

`list(mu = log-space mu, sigma = log-space sigma)`.

## See also

Other bayesian-fusion:
[`convert_beta_to_moments()`](https://jjmaynard.github.io/soilSIM/reference/convert_beta_to_moments.md),
[`convert_gamma_to_moments()`](https://jjmaynard.github.io/soilSIM/reference/convert_gamma_to_moments.md),
[`convert_lognormal_to_normal()`](https://jjmaynard.github.io/soilSIM/reference/convert_lognormal_to_normal.md),
[`convert_moments_to_beta()`](https://jjmaynard.github.io/soilSIM/reference/convert_moments_to_beta.md),
[`convert_moments_to_gamma()`](https://jjmaynard.github.io/soilSIM/reference/convert_moments_to_gamma.md),
[`fuse_beta()`](https://jjmaynard.github.io/soilSIM/reference/fuse_beta.md),
[`fuse_bivariate_normal()`](https://jjmaynard.github.io/soilSIM/reference/fuse_bivariate_normal.md),
[`fuse_distribution()`](https://jjmaynard.github.io/soilSIM/reference/fuse_distribution.md),
[`fuse_gamma()`](https://jjmaynard.github.io/soilSIM/reference/fuse_gamma.md),
[`fuse_normal_normal()`](https://jjmaynard.github.io/soilSIM/reference/fuse_normal_normal.md),
[`fuse_property()`](https://jjmaynard.github.io/soilSIM/reference/fuse_property.md),
[`update_prior()`](https://jjmaynard.github.io/soilSIM/reference/update_prior.md)
