# Determine candidate distributions to fit for a property

Computes the existing type-appropriate heuristic candidate list, then,
if `config$distribution_methods` is supplied, narrows to the
intersection of the two - keeping the heuristic as a soft prior rather
than silently discarding an explicit user request. If the intersection
is empty (the user asked only for families the heuristic wouldn't have
suggested), falls back to the user's list unfiltered, with a logged
warning, rather than ignoring it entirely.

## Usage

``` r
distributions_for_properties(property_name, values, config = NULL)
```

## Arguments

- property_name:

  Property name (drives the type-appropriate heuristic).

- values:

  Numeric vector of observed property values (currently unused by the
  heuristic itself, but taken for interface symmetry with callers and
  potential future data-driven refinement).

- config:

  Optional analysis configuration; reads `config$distribution_methods`.

## Value

Character vector of candidate distribution names.

## See also

Other statistics:
[`analyze_property_distributions()`](https://jjmaynard.github.io/soilSIM/reference/analyze_property_distributions.md),
[`analyze_soil_statistics()`](https://jjmaynard.github.io/soilSIM/reference/analyze_soil_statistics.md),
[`analyze_texture_correlations()`](https://jjmaynard.github.io/soilSIM/reference/analyze_texture_correlations.md),
[`compute_property_statistics()`](https://jjmaynard.github.io/soilSIM/reference/compute_property_statistics.md),
[`compute_stratified_correlations()`](https://jjmaynard.github.io/soilSIM/reference/compute_stratified_correlations.md),
[`default_statistics_config()`](https://jjmaynard.github.io/soilSIM/reference/default_statistics_config.md),
[`detect_comprehensive_outliers()`](https://jjmaynard.github.io/soilSIM/reference/detect_comprehensive_outliers.md),
[`diagnose_workflow()`](https://jjmaynard.github.io/soilSIM/reference/diagnose_workflow.md),
[`fit_property_distributions()`](https://jjmaynard.github.io/soilSIM/reference/fit_property_distributions.md),
[`generate_statistical_quality_report()`](https://jjmaynard.github.io/soilSIM/reference/generate_statistical_quality_report.md),
[`identify_numeric_soil_properties()`](https://jjmaynard.github.io/soilSIM/reference/identify_numeric_soil_properties.md),
[`run_comprehensive_correlation_analysis()`](https://jjmaynard.github.io/soilSIM/reference/run_comprehensive_correlation_analysis.md),
[`validate_statistical_config()`](https://jjmaynard.github.io/soilSIM/reference/validate_statistical_config.md),
[`validate_statistical_results()`](https://jjmaynard.github.io/soilSIM/reference/validate_statistical_results.md)
