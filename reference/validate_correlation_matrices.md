# Validate a set of correlation matrices

Real implementation using the shared
[`validate_correlation_matrix()`](https://jjmaynard.github.io/soilSIM/reference/validate_correlation_matrix.md)
(`distributions.R`) - previously a stub always returning `valid=TRUE`.

## Usage

``` r
validate_correlation_matrices(matrices, config)
```

## Arguments

- matrices:

  Named list keyed by method, each entry either a matrix or
  `list(matrix=, ...)` (matching
  [`run_comprehensive_correlation_analysis()`](https://jjmaynard.github.io/soilSIM/reference/run_comprehensive_correlation_analysis.md)'s
  `$matrices` shape).

- config:

  Unused; kept for interface compatibility with callers.

## Value

`list(valid=, errors=, warnings=)`.
