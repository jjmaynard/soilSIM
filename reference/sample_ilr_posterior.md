# Draw posterior composition samples from a fused/fitted ILR-space (mu, Sigma)

Guaranteed (by construction, via
[`ilr_inverse()`](https://jjmaynard.github.io/soilSIM/reference/ilr_inverse.md))
to sum to `total` and stay within `[0, total]` for every draw.

## Usage

``` r
sample_ilr_posterior(mu, Sigma, n = 1000, total = 100)
```

## Arguments

- mu:

  Length-2 mean vector.

- Sigma:

  2x2 covariance matrix.

- n:

  Number of samples to draw.

- total:

  Composition target sum.

## Value

An `n x 3` matrix (`clay`, `sand`, `silt`).
