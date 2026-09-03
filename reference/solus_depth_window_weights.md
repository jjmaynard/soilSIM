# Trapezoidal Depth-Average Weights over Native SOLUS Depth Slices

Returns the weights that turn SOLUS100's point predictions at the native
depth slices into an estimate of the **depth-average** of the property
over `[top_depth, bottom_depth]` - `sum(weights * slice_values)`
approximates `(1 / (b - a)) * integral_a^b f(d) dd`, where `f` is taken
piecewise-linear between native slices (with `rule = 2` constant
extrapolation outside the `[0, 150]` span SOLUS predicts over). This is
the raster analogue of the thickness-weighted window mean the SSURGO
prior side already computes
([`aggregate_depth_window_by_replicate()`](https://jjmaynard.github.io/soilSIM/reference/aggregate_depth_window_by_replicate.md)),
so the two sides of
[`fuse_property_adaptive()`](https://jjmaynard.github.io/soilSIM/reference/fuse_property_adaptive.md)
represent the same estimand.

## Usage

``` r
solus_depth_window_weights(
  top_depth,
  bottom_depth,
  available_slices = c(0, 5, 15, 30, 60, 100, 150)
)
```

## Arguments

- top_depth, bottom_depth:

  Numeric depth window bounds in cm (`bottom_depth > top_depth`).

- available_slices:

  Native SOLUS depth points (default `c(0, 5, 15, 30, 60, 100, 150)`).

## Value

A named numeric vector over `sort(unique(available_slices))` (names are
the slice depths as character), summing to 1. Zero for slices that do
not contribute.

## Details

Exact for a linear depth profile (reduces to the midpoint value,
matching
[`closest_solus_depth_slice()`](https://jjmaynard.github.io/soilSIM/reference/closest_solus_depth_slice.md));
unbiased to second order for curved profiles (OC decay, clay bulge, pH
trend). Cost is a single weighted raster-stack sum -
`length(available_slices)` is always small, so this scales with the
number of slices, not the number of pixels.

Windows lying entirely at/beyond the deepest slice (or above the
shallowest) fall back to the single nearest slice (constant
extrapolation).

## See also

[`closest_solus_depth_slice()`](https://jjmaynard.github.io/soilSIM/reference/closest_solus_depth_slice.md),
[`fetch_solus_low_pred_high()`](https://jjmaynard.github.io/soilSIM/reference/fetch_solus_low_pred_high.md)
