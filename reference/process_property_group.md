# Process Property Group Function

Per-group
([`dplyr::group_modify()`](https://dplyr.tidyverse.org/reference/group_map.html)-invoked)
wrapper for Strategies 1-3. Computes the problematic-cell mask for *this
group* - a group-modify reorders rows into single-group subsets, so a
whole-frame mask can't be passed in - then delegates to
[`infill_missing_property_data()`](https://jjmaynard.github.io/soilSIM/reference/infill_missing_property_data.md).

## Usage

``` r
process_property_group(
  group,
  property_name,
  property_config,
  max_depth,
  verbose
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

- verbose:

  Whether to provide progress messages

## Value

Data frame with infilled property data
