# Identify Soil Property Columns

Finds base soil-property names (a known SSURGO property with an `_r`
representative-value column present in `df`).

## Usage

``` r
identify_soil_property_columns_working(df)
```

## Arguments

- df:

  Horizon data frame.

## Value

Character vector of property base names (e.g. `"sandtotal"`).
