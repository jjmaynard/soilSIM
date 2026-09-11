# Property Data Cleaning

The single property-value cleaner for the whole package: parse
string/factor cells to numeric, null non-finite values, optionally flag
statistical outliers, then clamp to hardcoded physical-plausibility
bounds. Both the tabular pipeline
([`process_ssurgo_data()`](https://jjmaynard.github.io/soilSIM/reference/process_ssurgo_data.md))
and the infilling hierarchy
([`infill_soil_property()`](https://jjmaynard.github.io/soilSIM/reference/infill_soil_property.md))
call this.

## Usage

``` r
clean_property_data(
  df,
  property_name,
  outlier_policy = c("soil_aware", "aggressive_iqr", "none"),
  validation_config = NULL,
  generate_report = FALSE,
  verbose = FALSE
)
```

## Arguments

- df:

  Input data frame.

- property_name:

  Name of the property to clean.

- outlier_policy:

  How to flag statistical outliers on the `_r` column, before the
  physical-bound clamp:

  `"soil_aware"` (default)

  :   [`detect_statistical_outliers_soil_aware()`](https://jjmaynard.github.io/soilSIM/reference/detect_statistical_outliers_soil_aware.md) -
      flags only physically impossible values for pH / texture / bulk
      density, and a very conservative `IQR x 5` for everything else.
      The right default for a package whose job is uncertainty
      propagation: a real tail value is signal, not noise.

  `"aggressive_iqr"`

  :   generic `IQR x 3` for every property alike. Opt-in only.

  `"none"`

  :   skip outlier flagging; still parse, null non-finite, and clamp.

- validation_config:

  Optional per-property range rule (a `list(min=, max=)`-shaped entry);
  violations are counted into the report, values are not changed.

- generate_report:

  Whether to return a cleaning report alongside the data.

- verbose:

  Whether to provide progress messages.

## Value

`list(data = <cleaned df>)`, plus `report` when
`generate_report = TRUE`.

## See also

Other ssurgo-adapter:
[`add_range_rule()`](https://jjmaynard.github.io/soilSIM/reference/add_range_rule.md),
[`add_relationship_rule()`](https://jjmaynard.github.io/soilSIM/reference/add_relationship_rule.md),
[`apply_validation_rules()`](https://jjmaynard.github.io/soilSIM/reference/apply_validation_rules.md),
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
[`process_ssurgo_horizons()`](https://jjmaynard.github.io/soilSIM/reference/process_ssurgo_horizons.md),
[`property_contextual_ranges()`](https://jjmaynard.github.io/soilSIM/reference/property_contextual_ranges.md),
[`summarize_unsuitable_horizons()`](https://jjmaynard.github.io/soilSIM/reference/summarize_unsuitable_horizons.md)
