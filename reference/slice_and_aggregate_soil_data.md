# Slice and Aggregate Soil Data at Specified Depth Intervals

Slices a horizon data frame into 1 cm increments based on
`hzdept_r`/`hzdepb_r`, then averages every numeric column within each
requested depth range.

## Usage

``` r
slice_and_aggregate_soil_data(df, depth_ranges = list(c(0, 30), c(30, 100)))
```

## Arguments

- df:

  A data frame with `hzdept_r`/`hzdepb_r` depth-range columns.

- depth_ranges:

  A list of length-2 `c(top, bottom)` vectors.

## Value

A data frame with one row per depth range, mean values of soil
properties for each.
