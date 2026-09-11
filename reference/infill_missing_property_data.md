# Infill Missing Property Data

Applies the "within-group" recovery strategies, in order of reliability,
only using suitable horizons: horizon-name matching, depth-weighted
averaging, and within-component interpolation. This function is called
once per grouping unit (typically per `cokey`) via
[`process_property_group()`](https://jjmaynard.github.io/soilSIM/reference/process_property_group.md),
so it cannot see data from other groups - the group-spanning strategies
(cross-component interpolation, related-property estimation, and the
final group-mean fallback) are applied afterward at the whole-dataset
level by
[`infill_soil_property()`](https://jjmaynard.github.io/soilSIM/reference/infill_soil_property.md)
(see
[`cross_component_property_interpolation()`](https://jjmaynard.github.io/soilSIM/reference/cross_component_property_interpolation.md),
[`related_property_estimation()`](https://jjmaynard.github.io/soilSIM/reference/related_property_estimation.md),
and
[`apply_group_fallback_mean()`](https://jjmaynard.github.io/soilSIM/reference/apply_group_fallback_mean.md)).

## Usage

``` r
infill_missing_property_data(
  group,
  property_name,
  problematic_mask,
  property_config
)
```

## Arguments

- group:

  Data frame group to process

- property_name:

  Name of the property to infill

- problematic_mask:

  Logical vector for problematic horizons

- property_config:

  Property configuration

## Value

Data frame with infilled property data

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
[`infill_water_retention_saxton_rawls_integrated()`](https://jjmaynard.github.io/soilSIM/reference/infill_water_retention_saxton_rawls_integrated.md),
[`learn_property_ranges()`](https://jjmaynard.github.io/soilSIM/reference/learn_property_ranges.md),
[`process_soil_properties_comprehensive()`](https://jjmaynard.github.io/soilSIM/reference/process_soil_properties_comprehensive.md),
[`process_ssurgo_components()`](https://jjmaynard.github.io/soilSIM/reference/process_ssurgo_components.md),
[`process_ssurgo_data()`](https://jjmaynard.github.io/soilSIM/reference/process_ssurgo_data.md),
[`process_ssurgo_horizons()`](https://jjmaynard.github.io/soilSIM/reference/process_ssurgo_horizons.md),
[`property_contextual_ranges()`](https://jjmaynard.github.io/soilSIM/reference/property_contextual_ranges.md),
[`summarize_unsuitable_horizons()`](https://jjmaynard.github.io/soilSIM/reference/summarize_unsuitable_horizons.md)
