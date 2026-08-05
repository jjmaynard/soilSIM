# Repair a near-correlation matrix to be positive definite

Floors small/negative eigenvalues, reconstructs, and rescales back to a
unit-diagonal correlation matrix. Relocated verbatim from
`mod05_monte_carlo.R` (already correct there - not one of the confirmed
bugs, just moved so `statistics.R` can share it too).

## Usage

``` r
ensure_positive_definite_matrix(matrix, min_eigenvalue = 1e-06)
```

## Arguments

- matrix:

  A square numeric matrix.

- min_eigenvalue:

  Eigenvalue floor.

## Value

A positive-definite correlation matrix, same dimensions.
