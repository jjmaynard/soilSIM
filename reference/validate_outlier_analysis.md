# Validate an Enhanced Outlier Analysis Result

Flags any property/method whose `outlier_rate` exceeds a
config-defaulted threshold
(`config$statistical_analysis$outlier_rate_warn_threshold`, default
0.15 - a placeholder that should be tuned against real fixture data, not
treated as calibrated), and a multivariate step with a non-finite
`threshold` or `n_outliers`. Never errors.

## Usage

``` r
validate_outlier_analysis(outlier_analysis, config)
```

## Arguments

- outlier_analysis:

  Result of
  [`detect_comprehensive_outliers()`](https://jjmaynard.github.io/soilSIM/reference/detect_comprehensive_outliers.md)
  (or `_safe`).

- config:

  Analysis configuration.

## Value

`list(warnings = <character>)`.
