# Draw a Joint Depth x Property Gaussian Copula Sample

Draws `n_sims` realizations of an `n_depths x k` standard-normal field
whose ROWS (depths) are correlated according to `R_depth` and whose
COLUMNS (properties) are correlated according to `R_prop`,
SIMULTANEOUSLY - the standard separable/Kronecker-structured
multivariate-normal sampling identity
`Z = L_depth %*% Eps %*% t(L_prop)` (where `L_depth`/`L_prop` are
Cholesky factors of `R_depth`/`R_prop` and `Eps` is iid standard
normal), which achieves `Cov(vec(Z)) = R_prop (x) R_depth` (a Kronecker
product) without ever materializing that full
`(n_depths*k) x (n_depths*k)` matrix. Vectorized across all `n_sims`
realizations at once via array reshaping (no explicit per-realization
loop).

## Usage

``` r
sample_joint_depth_property_copula(R_depth, R_prop, n_sims, seed = NULL)
```

## Arguments

- R_depth:

  An `n_depths x n_depths` correlation matrix (e.g. from
  [`build_depth_correlation_kernel()`](https://jjmaynard.github.io/soilSIM/reference/build_depth_correlation_kernel.md)).

- R_prop:

  A `k x k` correlation matrix (the same flat property-correlation
  matrix
  [`simulate_correlated_triangular()`](https://jjmaynard.github.io/soilSIM/reference/simulate_correlated_triangular.md)
  already uses for within-horizon correlation).

- n_sims:

  Number of joint realizations to draw.

- seed:

  Optional integer seed for reproducibility.

## Value

A `c(n_depths, k, n_sims)` array of standard-normal values (mean 0,
variance 1 marginally), jointly correlated across both the depth and
property dimensions as specified.
