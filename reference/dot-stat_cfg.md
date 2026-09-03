# Resolve a statistical-analysis config value: nested location, then flat, then default

[`get_statistical_analysis_defaults()`](https://jjmaynard.github.io/soilSIM/reference/get_statistical_analysis_defaults.md)
nests every statistical parameter under
`config$statistical_analysis$...`, but the enhanced-chain functions
([`analyze_property_distributions()`](https://jjmaynard.github.io/soilSIM/reference/analyze_property_distributions.md),
[`detect_comprehensive_outliers()`](https://jjmaynard.github.io/soilSIM/reference/detect_comprehensive_outliers.md),
[`run_comprehensive_correlation_analysis()`](https://jjmaynard.github.io/soilSIM/reference/run_comprehensive_correlation_analysis.md))
were written reading them flat off `config$...` with no fallback - under
the package's own documented default usage
(`analyze_soil_statistics(horizon_data)`, no `analysis_config`) that
made `config$minimum_observations` `NULL` and hard-errored the enhanced
chain. This mirrors
[`analyze_soil_statistics()`](https://jjmaynard.github.io/soilSIM/reference/analyze_soil_statistics.md)'s
own `config$statistical_analysis$X %||% config$X %||% <default>` pattern
(see its `min_quality_score` line). Defaults here match the `_safe`
chain's hardcoded values and
[`get_statistical_analysis_defaults()`](https://jjmaynard.github.io/soilSIM/reference/get_statistical_analysis_defaults.md).

## Usage

``` r
.stat_cfg(config, key, default)
```

## Arguments

- config:

  The (merged) analysis configuration.

- key:

  Parameter name.

- default:

  Fallback when the key is absent in both locations.
