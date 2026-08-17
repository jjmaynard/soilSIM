# Safe Statistical Results Validation

Safe Statistical Results Validation

## Usage

``` r
validate_statistical_results_safe(
  correlation_analysis,
  distribution_analysis,
  outlier_analysis,
  property_statistics,
  config
)
```

## Arguments

- correlation_analysis:

  Correlation analysis results.

- distribution_analysis:

  Distribution analysis results.

- outlier_analysis:

  Outlier analysis results.

- property_statistics:

  Property statistics results.

- config:

  Analysis configuration.

## Value

`list(overall_valid=, errors=, warnings=, validation_score=, component_validations=)`.
