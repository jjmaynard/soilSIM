# Process Property Group Function

Process Property Group Function

## Usage

``` r
process_property_group(
  group,
  property_name,
  problematic_mask,
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

- problematic_mask:

  Logical vector for problematic horizons (IGNORED - recalculated)

- property_config:

  Property configuration

- max_depth:

  Maximum depth for processing

- verbose:

  Whether to provide progress messages

## Value

Data frame with infilled property data
