# Thickness-Weighted Mean of Simulated Properties over a Depth Window

For each replicate (a unique combination of `replicate_cols`), averages
every property column across that replicate's horizons, weighted by how
much of each horizon's thickness overlaps `[top_depth, bottom_depth]`.
Horizons with zero overlap contribute nothing.

## Usage

``` r
aggregate_depth_window_by_replicate(
  sim_long,
  top_depth,
  bottom_depth,
  property_cols,
  replicate_cols = c("mukey", "cokey", "simulation_number")
)
```

## Arguments

- sim_long:

  Long-format simulated data with `hzdept_r`/`hzdepb_r` and
  `replicate_cols`.

- top_depth, bottom_depth:

  Numeric depth window bounds in cm.

- property_cols:

  Character vector of property columns to aggregate.

- replicate_cols:

  Columns identifying one replicate (default
  `c("mukey", "cokey", "simulation_number")`).

## Value

One row per replicate, with `top`/`bottom` set to the requested window
and each property column replaced by its overlap-weighted mean (`NA` if
the replicate has no overlap with the window).
