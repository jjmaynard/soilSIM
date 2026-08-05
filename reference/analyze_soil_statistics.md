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
