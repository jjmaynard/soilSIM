# Snap a Depth Window to the Nearest Native SOLUS Depth Slice

[`soilDB::fetchSOLUS()`](http://ncss-tech.github.io/soilDB/reference/fetchSOLUS.md)'s
`depth_slices` are single depth points (default
`c(0, 5, 15, 30, 60, 100, 150)`, live-verified), not ranges, so a
`[top_depth, bottom_depth]` window must be approximated by one native
slice. This function did not exist anywhere in the source this was
ported from (referenced but never defined) - genuinely new design, not a
port: snaps to the native slice closest to the window's midpoint, the
simplest representative-point choice for an otherwise point-based data
source.

## Usage

``` r
closest_solus_depth_slice(
  top_depth,
  bottom_depth,
  available_slices = c(0, 5, 15, 30, 60, 100, 150)
)
```

## Arguments

- top_depth, bottom_depth:

  Numeric depth window bounds in cm.

- available_slices:

  Native SOLUS depth points to snap to.

## Value

A single numeric value from `available_slices`.
