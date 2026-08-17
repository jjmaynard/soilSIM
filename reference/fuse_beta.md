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
