# Depth-Trend GP Adjustment, Guarded by `GPfit` Availability

Applies
[`apply_local_gp_adjustments()`](https://jjmaynard.github.io/soilSIM/reference/apply_local_gp_adjustments.md)
(`R/multivariate-adjustment.R` - fits its own local GP per property from
each cokey's own within-simulation depth trend, no pre-supplied GP
models needed) per cokey, when `GPfit` is installed and a cokey has
enough distinct depths. Cokeys with fewer than `min_depths` distinct
depths pass through unadjusted, exactly as the original per-cokey guard
did.

## Usage

``` r
maybe_adjust_soil_data_depth_trend(sim_long, properties, min_depths = 2)
```

## Arguments

- sim_long:

  Long-format simulated data with `cokey`, `hzdept_r`, and property
  columns.

- properties:

  Character vector of property column names to adjust.

- min_depths:

  Minimum distinct depths required to attempt GP fitting (default 2,
  matching the source's `length(unique_depths) >= 2` guard).

## Value

`sim_long`, depth-trend-adjusted where possible.
