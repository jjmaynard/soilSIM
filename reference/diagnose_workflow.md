# Validate Complete Workflow

Assesses the entire soil simulation workflow end to end.

## Usage

``` r
diagnose_workflow(
  workflow_results,
  original_data = NULL,
  validation_config = NULL,
  generate_plots = TRUE,
  output_dir = NULL,
  verbose = getOption("ssurgo.verbose", FALSE)
)
```

## Arguments

- workflow_results:

  Complete workflow results including all intermediate steps

- original_data:

  Original SSURGO/NRCS data for comparison

- validation_config:

  Configuration for validation parameters (uses package defaults)

- generate_plots:

  Whether to generate diagnostic plots (default = TRUE)

- output_dir:

  Directory for saving validation outputs

- verbose:

  Logical; if `TRUE`, temporarily raises the package's log level so
  `INFO`-level progress messages print for the duration of this call
  (default `FALSE` - quiet). See
  [`set_verbose_logging()`](https://jjmaynard.github.io/soilSIM/reference/set_verbose_logging.md).

## Value

Comprehensive validation results, as a `soilSIM_diagnostics` object

## See also

Other statistics:
[`analyze_property_distributions()`](https://jjmaynard.github.io/soilSIM/reference/analyze_property_distributions.md),
[`analyze_soil_statistics()`](https://jjmaynard.github.io/soilSIM/reference/analyze_soil_statistics.md),
[`analyze_texture_correlations()`](https://jjmaynard.github.io/soilSIM/reference/analyze_texture_correlations.md),
[`compute_property_statistics()`](https://jjmaynard.github.io/soilSIM/reference/compute_property_statistics.md),
[`compute_stratified_correlations()`](https://jjmaynard.github.io/soilSIM/reference/compute_stratified_correlations.md),
[`default_statistics_config()`](https://jjmaynard.github.io/soilSIM/reference/default_statistics_config.md),
[`detect_comprehensive_outliers()`](https://jjmaynard.github.io/soilSIM/reference/detect_comprehensive_outliers.md),
[`distributions_for_properties()`](https://jjmaynard.github.io/soilSIM/reference/distributions_for_properties.md),
[`fit_property_distributions()`](https://jjmaynard.github.io/soilSIM/reference/fit_property_distributions.md),
[`generate_statistical_quality_report()`](https://jjmaynard.github.io/soilSIM/reference/generate_statistical_quality_report.md),
[`identify_numeric_soil_properties()`](https://jjmaynard.github.io/soilSIM/reference/identify_numeric_soil_properties.md),
[`run_comprehensive_correlation_analysis()`](https://jjmaynard.github.io/soilSIM/reference/run_comprehensive_correlation_analysis.md),
[`validate_statistical_config()`](https://jjmaynard.github.io/soilSIM/reference/validate_statistical_config.md),
[`validate_statistical_results()`](https://jjmaynard.github.io/soilSIM/reference/validate_statistical_results.md)
