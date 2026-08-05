# Fuse two independent bivariate Normal beliefs about the same 2D quantity

Precision MATRICES add - the direct multivariate generalization of
[`bayes_update_normal_normal()`](https://jjmaynard.github.io/soilSIM/reference/bayes_update_normal_normal.md)'s
scalar-precision addition.

## Usage

``` r
fuse_bivariate_normal(mu1, Sigma1, mu2, Sigma2)
```

## Arguments

- mu1, mu2:

  Length-2 mean vectors.

- Sigma1, Sigma2:

  2x2 covariance matrices.

## Value

`list(mu = fused mean vector, Sigma = fused covariance matrix)`.
