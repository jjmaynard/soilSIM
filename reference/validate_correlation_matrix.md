# Validate a correlation matrix's shape and positive-definiteness

Relocated verbatim from `mod05_monte_carlo.R`.

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
