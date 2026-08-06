# Depth-Trend GP Adjustment for a Single Cokey's Data

The per-cokey unit of work
[`maybe_adjust_soil_data_depth_trend()`](https://jjmaynard.github.io/soilSIM/reference/maybe_adjust_soil_data_depth_trend.md)
maps over every cokey, factored out so it can be dispatched to parallel
workers unchanged.

## Usage

``` r
adjust_one_cokey_depth_trend(cokey_data, properties, min_depths)
```

## Arguments

- cokey_data:

  Simulation data for a single cokey.

- properties:

  Character vector of property column names to adjust.

- min_depths:

  Minimum distinct depths required to attempt GP fitting.

## Value

`cokey_data`, depth-trend-adjusted where possible.
