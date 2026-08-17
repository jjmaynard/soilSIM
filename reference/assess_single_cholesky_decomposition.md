# Assess a Single Cholesky Decomposition

Actually attempts a Cholesky decomposition of `matrix_data` and measures
reconstruction error (`||L'L - matrix_data||` in Frobenius norm).

## Usage

``` r
assess_single_cholesky_decomposition(matrix_data, criteria)
```

## Arguments

- matrix_data:

  Matrix to decompose.

- criteria:

  List, optionally with `reconstruction_tolerance` (default 1e-10).

## Value

List with `decomposition_successful`, `reconstruction_error`,
`numerical_stability`.
