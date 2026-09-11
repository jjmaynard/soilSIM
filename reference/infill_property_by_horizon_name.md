# Horizon Name Property Infill

Infills missing values using horizons with matching or similar names,
only using suitable horizons as sources.

## Usage

``` r
infill_property_by_horizon_name(group, property_col, problematic_mask)
```

## Arguments

- group:

  Data frame group

- property_col:

  Property column name

- problematic_mask:

  Logical mask for problematic horizons

## Value

Data frame with infilled values
