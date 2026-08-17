# Simulate Correlated Samples from Triangular Distributions

Generates `n` correlated samples from `length(params)` triangular
distributions via a Cholesky decomposition of `correlation_matrix`
applied to standard-normal draws, then transforms each column to its own
triangular distribution via the inverse CDF.

## Usage

``` r
simulate_correlated_triangular(
  n,
  params,
  correlation_matrix,
  random_seed = NULL
)
```

## Arguments

- n:

  Integer, number of samples to generate.

- params:

  List of `c(a, b, c)` triples, one per distribution - **note the order
  here is (lower, mode, upper)**, unlike
  [`tri_dist()`](https://jjmaynard.github.io/soilSIM/reference/tri_dist.md)'s
  own `(a = lower, b = upper, c = mode)` convention. Preserved exactly
  as the legacy source defines it; a real source of confusion if assumed
  to match
  [`tri_dist()`](https://jjmaynard.github.io/soilSIM/reference/tri_dist.md)'s
  argument order.

- correlation_matrix:

  A square, positive-semi-definite correlation matrix, `length(params)`
  x `length(params)`.

- random_seed:

  Optional integer seed for reproducibility.

## Value

A matrix of correlated samples, `n` rows x `length(params)` columns.
