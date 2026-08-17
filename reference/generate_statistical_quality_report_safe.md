# Safe Quality Report Generation

Safe Quality Report Generation

## Usage

``` r
generate_statistical_quality_report_safe(
  original_data,
  processed_data,
  correlation_analysis,
  distribution_analysis,
  outlier_analysis,
  property_statistics,
  validation_results,
  data_validation,
  config
)
```

## Arguments

- original_data:

  Original input data.

- processed_data:

  Processed data.

- correlation_analysis:

  Correlation analysis results.

- distribution_analysis:

  Distribution analysis results.

- outlier_analysis:

  Outlier analysis results.

- property_statistics:

  Property statistics results.

- validation_results:

  Validation results.

- data_validation:

  Data validation results.

- config:

  Analysis configuration.

## Value

`list(overall_quality_score=, data_quality=, analysis_quality=, validation_summary=, recommendations=)`.
