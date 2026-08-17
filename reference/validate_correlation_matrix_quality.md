# Validate Correlation Matrix Quality

Inspects every matrix in `correlation_matrices` (a named list) for
positive-definiteness (via eigenvalues) and condition number, reusing
[`assess_correlation_matrix_properties()`](https://jjmaynard.github.io/soilSIM/reference/assess_correlation_matrix_properties.md)
per matrix and aggregating across all of them.

## Usage

``` r
validate_correlation_matrix_quality(correlation_matrices, criteria)
```

## Arguments

- correlation_matrices:

  Named list of correlation matrices.

- criteria:

  List, optionally with `condition_number_threshold` (default 1e12).

## Value

List with `positive_definite` (TRUE only if all matrices are),
`condition_number` (max across matrices), `matrix_quality`.
