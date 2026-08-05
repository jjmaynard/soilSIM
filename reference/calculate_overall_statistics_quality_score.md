# Calculate Overall Statistics Quality Score

Real composite score averaging whichever of `data_validation`'s quality
score, `validation_results`'s validation score, and the
[`assess_correlation_quality()`](https://jjmaynard.github.io/soilSIM/reference/assess_correlation_quality.md)/[`assess_distribution_quality()`](https://jjmaynard.github.io/soilSIM/reference/assess_distribution_quality.md)/
[`assess_outlier_quality()`](https://jjmaynard.github.io/soilSIM/reference/assess_outlier_quality.md)
scores are available.

## Usage

``` r
calculate_overall_statistics_quality_score(
  data_validation,
  correlation_analysis,
  distribution_analysis,
  outlier_analysis,
  validation_results
)
```

## Arguments

- data_validation:

  Result of
  [`validate_data_quality()`](https://jjmaynard.github.io/soilSIM/reference/validate_data_quality.md),
  with an `overall_quality$score` field.

- correlation_analysis:

  Result of correlation analysis (see
  [`assess_correlation_quality()`](https://jjmaynard.github.io/soilSIM/reference/assess_correlation_quality.md)),
  or `NULL`.

- distribution_analysis:

  Result of distribution analysis (see
  [`assess_distribution_quality()`](https://jjmaynard.github.io/soilSIM/reference/assess_distribution_quality.md)),
  or `NULL`.

- outlier_analysis:

  Result of
  [`detect_comprehensive_outliers()`](https://jjmaynard.github.io/soilSIM/reference/detect_comprehensive_outliers.md),
  or `NULL`.

- validation_results:

  Result of statistical-results validation, with a `validation_score`
  field.

## Value

Single numeric quality score in `[0, 1]`, or `NA` if nothing was
available.
