# Process Property Group Fallback Function

Companion to
[`process_property_group()`](https://jjmaynard.github.io/soilSIM/reference/process_property_group.md):
applies
[`apply_group_fallback_mean()`](https://jjmaynard.github.io/soilSIM/reference/apply_group_fallback_mean.md)
(Strategy 6, last resort) per grouping unit, after the whole-dataset
Strategy 4/5 passes (cross-component interpolation, related-property
estimation) have already run in
[`infill_soil_property()`](https://jjmaynard.github.io/soilSIM/reference/infill_soil_property.md).
Recalculates its own problematic mask the same way
[`process_property_group()`](https://jjmaynard.github.io/soilSIM/reference/process_property_group.md)
does, so it only fills cells that are still genuinely missing at this
point.

## Usage

``` r
process_property_group_fallback(
  group,
  property_name,
  property_config,
  max_depth
)
```

## Arguments

- group:

  Data frame group to process

- property_name:

  Name of the property to infill

- property_config:

  Property configuration

- max_depth:

  Maximum depth for processing

## Value

Data frame with fallback-filled property data
