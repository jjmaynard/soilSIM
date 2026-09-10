# Validate a correlation matrix's shape and positive-definiteness

Validate a correlation matrix's shape and positive-definiteness

## Usage

``` r
validate_correlation_matrix(corr_matrix, properties)
```

## Arguments

- corr_matrix:

  A candidate correlation matrix.

- properties:

  Character vector the matrix's dimensions should match.

## Value

`list(valid=, message=)`.
