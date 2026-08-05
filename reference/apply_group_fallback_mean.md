# Apply Group-Mean Fallback

Strategy 6 (last resort) of the property recovery hierarchy: fills any
still-missing suitable-horizon values with a depth-weighted (or plain)
mean computed from the suitable horizons of `group`. Extracted as its
own function so it can be applied after the group-spanning Strategy 4/5
passes (cross-component interpolation, related-property estimation) have
already had a chance to fill values - see
[`infill_soil_property()`](https://jjmaynard.github.io/soilSIM/reference/infill_soil_property.md).

## Usage

``` r
apply_group_fallback_mean(group, property_col, problematic_mask)
```

## Arguments

- group:

  Data frame group to process

- property_col:

  Name of the property column to infill

- problematic_mask:

  Logical vector for problematic horizons

## Value

Data frame with fallback-filled values
