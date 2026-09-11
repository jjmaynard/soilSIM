# Summarize Unsuitable Horizons

Logs and returns a summary of which horizon types were excluded from
infilling (the flagging itself is handled by
[`is_unsuitable()`](https://jjmaynard.github.io/soilSIM/reference/is_unsuitable.md)).

## Usage

``` r
summarize_unsuitable_horizons(df, hzname_col = "hzname")
```

## Arguments

- df:

  Data frame with an `unsuitable_horizon` logical column (see
  [`is_unsuitable()`](https://jjmaynard.github.io/soilSIM/reference/is_unsuitable.md))

- hzname_col:

  Name of the horizon-name column

## Value

List with `n_unsuitable` (count) and `horizon_types` (unique excluded
horizon names)

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
[`process_ssurgo_horizons()`](https://jjmaynard.github.io/soilSIM/reference/process_ssurgo_horizons.md),
[`property_contextual_ranges()`](https://jjmaynard.github.io/soilSIM/reference/property_contextual_ranges.md)
