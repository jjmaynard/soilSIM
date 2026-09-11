# Process Horizon Data (Working Compatible)

Process Horizon Data (Working Compatible)

## Usage

``` r
process_ssurgo_horizons(
  raw_data,
  detect_unsuitable = TRUE,
  advanced_cleaning = TRUE,
  standardize_names = TRUE,
  remove_invalid = TRUE,
  calculate_derived = TRUE,
  max_depth = DEFAULT_MAX_DEPTH_CM,
  verbose = FALSE
)
```

## Arguments

- raw_data:

  Raw SSURGO data containing horizon information

- detect_unsuitable:

  Logical; detect unsuitable horizons

- advanced_cleaning:

  Logical; apply advanced data cleaning functions

- standardize_names:

  Logical; standardize column names

- remove_invalid:

  Logical; remove invalid horizon records

- calculate_derived:

  Logical; calculate derived properties

- max_depth:

  Maximum depth for processing

- verbose:

  Logical; provide progress messages

## Value

List with processed horizon data and processing statistics

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
[`process_ssurgo_data()`](https://jjmaynard.github.io/soilSIM/reference/process_ssurgo_data.md),
[`property_contextual_ranges()`](https://jjmaynard.github.io/soilSIM/reference/property_contextual_ranges.md),
[`summarize_unsuitable_horizons()`](https://jjmaynard.github.io/soilSIM/reference/summarize_unsuitable_horizons.md)
