# Cross-Component Property Interpolation

Fills missing values in suitable horizons using data from other soil
components (other `cokey` groups within the same `group`) at similar
depths, weighted by depth proximity. Only suitable horizons are used as
sources.

## Usage

``` r
cross_component_property_interpolation(group, property_col)
```

## Arguments

- group:

  Data frame group to process

- property_col:

  Name of the property column to infill

## Value

Data frame with infilled values
