# Assess Correlation Matrix Properties

Checks a correlation matrix's positive-definiteness (via eigenvalues)
and condition number.

## Usage

``` r
assess_correlation_matrix_properties(cor_matrix)
```

## Arguments

- cor_matrix:

  Correlation matrix.

## Value

List with `is_positive_definite`, `condition_number`,
`eigenvalue_range`, `matrix_stability`.
