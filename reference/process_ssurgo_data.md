# Process SSURGO Data (Main Entry Point)

Version that maintains compatibility with working infill_soil_property()
while adding comprehensive processing, validation, and reporting
capabilities. Uses shared utilities for general processing tasks.

## Usage

``` r
process_ssurgo_data(
  raw_data,
  processing_options = list(),
  validate_results = TRUE,
  max_depth = DEFAULT_MAX_DEPTH_CM,
  verbose = getOption("ssurgo.verbose", FALSE)
)
```

## Arguments

- raw_data:

  Raw SSURGO data from download functions

- processing_options:

  List of processing options and parameters

- validate_results:

  Logical; validate processing results (default: TRUE)

- max_depth:

  Maximum depth (cm) for processing (default: 250)

- verbose:

  Logical; provide detailed progress messages

## Value

List containing:

- processed_data: Combined processed data (compatible with infill
  functions)

- horizon_data: Processed horizon data

- component_data: Processed component data

- processing_metadata: Processing information and statistics

- validation_results: Data validation results

- quality_report: Comprehensive quality assessment

## See also

Other ssurgo-adapter:
[`add_range_rule()`](https://jjmaynard.github.io/soilSIM/reference/add_range_rule.md),
[`add_relationship_rule()`](https://jjmaynard.github.io/soilSIM/reference/add_relationship_rule.md),
[`apply_validation_rules()`](https://jjmaynard.github.io/soilSIM/reference/apply_validation_rules.md),
[`clean_property_data()`](https://jjmaynard.github.io/soilSIM/reference/clean_property_data.md),
[`create_custom_property_config()`](https://jjmaynard.github.io/soilSIM/reference/create_custom_property_config.md),
[`create_validation_config()`](https://jjmaynard.github.io/soilSIM/reference/create_validation_config.md),
[`default_property_config()`](https://jjmaynard.github.io/soilSIM/reference/default_property_config.md),
[`impute_rfv_values()`](https://jjmaynard.github.io/soilSIM/reference/impute_rfv_values.md),
[`infill_missing_property_data()`](https://jjmaynard.github.io/soilSIM/reference/infill_missing_property_data.md),
[`infill_water_retention_saxton_rawls_integrated()`](https://jjmaynard.github.io/soilSIM/reference/infill_water_retention_saxton_rawls_integrated.md),
[`learn_property_ranges()`](https://jjmaynard.github.io/soilSIM/reference/learn_property_ranges.md),
[`process_soil_properties_comprehensive()`](https://jjmaynard.github.io/soilSIM/reference/process_soil_properties_comprehensive.md),
[`process_ssurgo_components()`](https://jjmaynard.github.io/soilSIM/reference/process_ssurgo_components.md),
[`process_ssurgo_horizons()`](https://jjmaynard.github.io/soilSIM/reference/process_ssurgo_horizons.md),
[`property_contextual_ranges()`](https://jjmaynard.github.io/soilSIM/reference/property_contextual_ranges.md),
[`summarize_unsuitable_horizons()`](https://jjmaynard.github.io/soilSIM/reference/summarize_unsuitable_horizons.md)
