# Snap a Depth Window to the Nearest Native SOLUS Depth Slice

[`soilDB::fetchSOLUS()`](http://ncss-tech.github.io/soilDB/reference/fetchSOLUS.md)'s
`depth_slices` are single depth points (default
`c(0, 5, 15, 30, 60, 100, 150)`, live-verified), not ranges, so a
`[top_depth, bottom_depth]` window must be approximated. This function
did not exist anywhere in the source this was ported from (referenced
but never defined) - genuinely new design, not a port: snaps to the
native slice closest to the window's midpoint.

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

## Superseded for the fusion pipeline

[`fetch_solus_low_pred_high()`](https://jjmaynard.github.io/soilSIM/reference/fetch_solus_low_pred_high.md)
no longer uses this - a single midpoint point-value only equals the
window depth-average when the property is linear in depth, which most
soil properties are not, and it left the SOLUS likelihood representing a
different quantity (a point) than the SSURGO prior side (a genuine
thickness-weighted window mean, see
[`aggregate_depth_window_by_replicate()`](https://jjmaynard.github.io/soilSIM/reference/aggregate_depth_window_by_replicate.md)).
It now uses
[`solus_depth_window_weights()`](https://jjmaynard.github.io/soilSIM/reference/solus_depth_window_weights.md)
for a trapezoidal depth-average across the native slices spanning the
window. This function is retained as public API for back-compat and for
callers that genuinely want a single nearest slice.

## See also

[`solus_depth_window_weights()`](https://jjmaynard.github.io/soilSIM/reference/solus_depth_window_weights.md)
