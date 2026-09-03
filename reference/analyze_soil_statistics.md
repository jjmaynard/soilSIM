# Comprehensive Statistical Analysis (Main Entry Point)

Enhanced statistical analysis leveraging Module 0 utilities for data
validation, error handling, logging, and statistical computations.

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

  Processed SSURGO data from Module 2

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
  flip is a separate, later, benchmark-driven decision (see
  `RASTER_STATISTICS_INTEGRATION_PLAN.md` Problem B).

- verbose:

  Logical; provide detailed progress messages

## Value

List containing comprehensive statistical analysis results

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
