# Analyze texture (sand/silt/clay) correlations, in both raw and ILR space

Raw Pearson correlations among compositional (simplex-constrained) parts
like sand/silt/clay percentages are spuriously negative due to the
sum-to-100 constraint - `ilr_correlations` (isometric log-ratio space,
via the dependency-free
[`ilr_forward()`](https://jjmaynard.github.io/soilSIM/reference/ilr_forward.md)
in `core-distributions.R`) avoids that artifact and is the statistically
defensible view for compositional data. Both are reported:
`raw_correlations` for continuity with prior behavior,
`ilr_correlations` as the added analytical value (only computed when all
three fractions are present with at least 5 complete observations).

## Usage

``` r
analyze_texture_correlations(data, texture_properties, methods, config)
```

## Arguments

- data:

  Input data with texture `_r` columns.

- texture_properties:

  Character vector of available texture column names (e.g. a subset/all
  of `c("sandtotal_r","silttotal_r","claytotal_r")`).

- methods:

  Correlation methods (e.g. `c("pearson","spearman")`).

- config:

  Unused; kept for interface compatibility with callers.

## Value

`list(raw_correlations=, ilr_correlations=, n_observations=, note=)`.

## See also

Other statistics:
[`analyze_property_distributions()`](https://jjmaynard.github.io/soilSIM/reference/analyze_property_distributions.md),
[`analyze_soil_statistics()`](https://jjmaynard.github.io/soilSIM/reference/analyze_soil_statistics.md),
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
