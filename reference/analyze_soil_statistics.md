# Comprehensive Statistical Analysis (Main Entry Point)

Statistical analysis for data validation, error handling, logging, and
statistical computations.

## Usage

``` r
analyze_soil_statistics(
  processed_data,
  analysis_config = list(),
  correlation_methods = c("pearson", "spearman"),
  distribution_fitting = TRUE,
  outlier_detection = TRUE,
  validate_results = TRUE,
  use_enhanced_chain = FALSE,
  verbose = getOption("ssurgo.verbose", FALSE)
)
```

## Arguments

- processed_data:

  Processed SSURGO data from the statistics step

- analysis_config:

  List of analysis configuration parameters

- correlation_methods:

  Character vector of correlation methods to apply

- distribution_fitting:

  Logical; perform distribution fitting analysis

- outlier_detection:

  Logical; perform statistical outlier detection

- validate_results:

  Logical; validate statistical results (default: TRUE)

- use_enhanced_chain:

  Logical (default `FALSE`). `FALSE` runs the `_safe` chain -
  [`run_comprehensive_correlation_analysis_safe()`](https://jjmaynard.github.io/soilSIM/reference/run_comprehensive_correlation_analysis_safe.md)
  /
  [`analyze_property_distributions_safe()`](https://jjmaynard.github.io/soilSIM/reference/analyze_property_distributions_safe.md)
  /
  [`detect_comprehensive_outliers_safe()`](https://jjmaynard.github.io/soilSIM/reference/detect_comprehensive_outliers_safe.md)
  /
  [`validate_statistical_results_safe()`](https://jjmaynard.github.io/soilSIM/reference/validate_statistical_results_safe.md)
  /
  [`generate_statistical_quality_report_safe()`](https://jjmaynard.github.io/soilSIM/reference/generate_statistical_quality_report_safe.md) -
  preserving today's output shape exactly (pinned by an
  [`identical()`](https://rdrr.io/r/base/identical.html) regression
  test). `TRUE` routes through the richer enhanced functions instead
  (parametric distribution fits with AIC/BIC/convergence, per-method
  outlier detection, multivariate Mahalanobis outliers, stratified + ILR
  texture correlations). The two return **different, non-interchangeable
  shapes** - notably
  `distribution_analysis$fitted_distributions[[prop]]` is a flat stats
  list under `_safe` and a family-name-keyed list of fit objects under
  enhanced. Both chains stay public API. Whether the default should ever
  flip is a separate, later, benchmark-driven decision.

- verbose:

  Logical; provide detailed progress messages

## Value

List containing comprehensive statistical analysis results

## See also

Other statistics:
[`analyze_property_distributions()`](https://jjmaynard.github.io/soilSIM/reference/analyze_property_distributions.md),
[`analyze_texture_correlations()`](https://jjmaynard.github.io/soilSIM/reference/analyze_texture_correlations.md),
[`compute_property_statistics()`](https://jjmaynard.github.io/soilSIM/reference/compute_property_statistics.md),
[`compute_stratified_correlations()`](https://jjmaynard.github.io/soilSIM/reference/compute_stratified_correlations.md),
[`default_statistics_config()`](https://jjmaynard.github.io/soilSIM/reference/default_statistics_config.md),
[`detect_comprehensive_outliers()`](https://jjmaynard.github.io/soilSIM/reference/detect_comprehensive_outliers.md),
[`diagnose_workflow()`](https://jjmaynard.github.io/soilSIM/reference/diagnose_workflow.md),
[`distributions_for_properties()`](https://jjmaynard.github.io/soilSIM/reference/distributions_for_properties.md),
[`fit_property_distributions()`](https://jjmaynard.github.io/soilSIM/reference/fit_property_distributions.md),
[`generate_statistical_quality_report()`](https://jjmaynard.github.io/soilSIM/reference/generate_statistical_quality_report.md),
[`identify_numeric_soil_properties()`](https://jjmaynard.github.io/soilSIM/reference/identify_numeric_soil_properties.md),
[`run_comprehensive_correlation_analysis()`](https://jjmaynard.github.io/soilSIM/reference/run_comprehensive_correlation_analysis.md),
[`validate_statistical_config()`](https://jjmaynard.github.io/soilSIM/reference/validate_statistical_config.md),
[`validate_statistical_results()`](https://jjmaynard.github.io/soilSIM/reference/validate_statistical_results.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic statistical analysis
stats_result <- analyze_soil_statistics(
  processed_data,
  analysis_config = list(
    stratify_by_horizon = TRUE,
    include_texture_analysis = TRUE,
    correlation_threshold = 0.05
  )
)
} # }
```
