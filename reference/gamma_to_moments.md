# Gamma's mean/variance as a function of its own (shape, rate)

The inverse direction of
[`moments_to_gamma()`](https://jjmaynard.github.io/soilSIM/reference/moments_to_gamma.md) -
used to re-express an infeasible same-family fusion's inputs as Normal
moments before falling back to
[`bayes_update_normal_normal()`](https://jjmaynard.github.io/soilSIM/reference/bayes_update_normal_normal.md).

## Usage

``` r
gamma_to_moments(shape, rate)
```

## Arguments

- shape, rate:

  Gamma shape/rate parameters.

## Value

`list(mean=, var=)`.
