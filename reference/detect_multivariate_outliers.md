# Detect Multivariate Outliers

Real Mahalanobis-distance-based multivariate outlier detection across
`properties`, guarding a non-positive-definite covariance matrix with
[`Matrix::nearPD()`](https://rdrr.io/pkg/Matrix/man/nearPD.html)
(already an Import, no new dependency).

## Usage

``` r
detect_multivariate_outliers(data, properties, config)
```

## Arguments

- data:

  Input data frame.

- properties:

  Character vector of column names to include.

- config:

  List, optionally with `mahalanobis_alpha` (default 0.975, i.e. flag
  the outer 2.5% via a chi-squared cutoff).

## Value

List with `outliers` (logical vector, `NA` for incomplete rows),
`n_outliers`, `method`, `threshold`.
