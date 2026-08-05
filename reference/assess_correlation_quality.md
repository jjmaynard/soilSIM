# Assess Correlation Analysis Quality

Real score based on the fraction of correlation matrices that produced a
validation warning (see
[`validate_correlation_matrices()`](https://jjmaynard.github.io/soilSIM/reference/validate_correlation_matrices.md)'s
`$matrices`/`$validation$warnings` shape).

## Usage

``` r
assess_correlation_quality(correlation_analysis)
```

## Arguments

- correlation_analysis:

  List with `matrices` and `validation$warnings`.

## Value

List with `quality_score`, `n_matrices`, `n_warnings`.
