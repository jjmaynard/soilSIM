# Snap a Depth Window to the Nearest Native SOLUS Depth Slice

[`soilDB::fetchSOLUS()`](http://ncss-tech.github.io/soilDB/reference/fetchSOLUS.md)'s
`depth_slices` are single depth points (default
`c(0, 5, 15, 30, 60, 100, 150)`, live-verified), not ranges, so a
`[top_depth, bottom_depth]` window must be approximated. This function
snaps to the native slice closest to the window's midpoint.

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

## Relationship to the fusion pipeline

[`fetch_solus_low_pred_high()`](https://jjmaynard.github.io/soilSIM/reference/fetch_solus_low_pred_high.md)
uses
[`solus_depth_window_weights()`](https://jjmaynard.github.io/soilSIM/reference/solus_depth_window_weights.md)
for a trapezoidal depth-average across the native slices spanning the
window, so the SOLUS likelihood and the SSURGO prior both represent a
thickness-weighted window mean. This function is public API for callers
that want a single nearest slice.

## See also

[`solus_depth_window_weights()`](https://jjmaynard.github.io/soilSIM/reference/solus_depth_window_weights.md)
