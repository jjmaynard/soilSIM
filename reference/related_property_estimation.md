# Related Property Estimation

Estimates missing values from pedologically-related properties, using
property-type-specific relationships (e.g. texture sum-to-100
constraint, clay/organic-matter based CEC estimation, clay-based water
retention, horizon/organic-matter adjusted pH, depth/horizon/clay
adjusted organic matter, texture-adjusted bulk density). Only suitable
horizons are estimated.

## Usage

``` r
related_property_estimation(group, property_name, property_config)
```

## Arguments

- group:

  Data frame group to process

- property_name:

  Name of the property to estimate

- property_config:

  Property configuration list

## Value

Data frame with estimated values
