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
