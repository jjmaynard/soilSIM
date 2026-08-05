# Assess Distributional Coverage

Assesses `values`' own distributional spread: outlier percentage (IQR
method) and whether the interquartile range is non-degenerate.

## Usage

``` r
assess_distributional_coverage(values, criteria)
```

## Arguments

- values:

  Values to assess.

- criteria:

  List, optionally with `max_outlier_percentage` (default 10).

## Value

List with `distribution_coverage`, `outlier_percentage`,
`coverage_adequate`.
