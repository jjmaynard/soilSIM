# Assess Outlier Analysis Quality

Real score based on the overall outlier rate implied by
[`detect_comprehensive_outliers()`](https://jjmaynard.github.io/soilSIM/reference/detect_comprehensive_outliers.md)'s
`summary`/`property_outliers` shape.

## Usage

``` r
assess_outlier_quality(outlier_analysis)
```

## Arguments

- outlier_analysis:

  Result of
  [`detect_comprehensive_outliers()`](https://jjmaynard.github.io/soilSIM/reference/detect_comprehensive_outliers.md).

## Value

List with `quality_score`, `total_outliers`, `outlier_rate`.
