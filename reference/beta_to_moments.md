# Beta's mean/variance as a function of its own (alpha, beta)

The inverse direction of
[`moments_to_beta()`](https://jjmaynard.github.io/soilSIM/reference/moments_to_beta.md) -
used to re-express an infeasible same-family fusion's inputs as Normal
moments before falling back to
[`bayes_update_normal_normal()`](https://jjmaynard.github.io/soilSIM/reference/bayes_update_normal_normal.md).

## Usage

``` r
beta_to_moments(alpha, beta)
```

## Arguments

- alpha, beta:

  Beta shape parameters.

## Value

`list(mean=, var=)`.
