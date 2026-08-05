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
