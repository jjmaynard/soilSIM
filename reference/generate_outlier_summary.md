# Generate Outlier Summary

Aggregates
[`detect_comprehensive_outliers()`](https://jjmaynard.github.io/soilSIM/reference/detect_comprehensive_outliers.md)'s
`property_outliers` (per property, per method) and
`multivariate_outliers` into overall counts.

## Usage

``` r
generate_outlier_summary(outlier_results)
```

## Arguments

- outlier_results:

  Result of
  [`detect_comprehensive_outliers()`](https://jjmaynard.github.io/soilSIM/reference/detect_comprehensive_outliers.md)
  (or a list with the same `property_outliers`/`multivariate_outliers`
  shape).

## Value

List with `total_outliers`, `property_outlier_counts`,
`multivariate_outlier_count`.
