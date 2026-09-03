# Fetch SOLUS100 Low/Prediction/High Rasters for One Variable and Depth Window

Fetches every native SOLUS depth slice spanning
`[top_depth, bottom_depth]` and reduces them to a single
**depth-average** raster per output type via
[`solus_depth_window_weights()`](https://jjmaynard.github.io/soilSIM/reference/solus_depth_window_weights.md) -
a trapezoidal integral of SOLUS's point predictions over the window. See
[`solus_depth_window_weights()`](https://jjmaynard.github.io/soilSIM/reference/solus_depth_window_weights.md)
for why this replaces the previous single-nearest-slice snap
([`closest_solus_depth_slice()`](https://jjmaynard.github.io/soilSIM/reference/closest_solus_depth_slice.md)):
a point value only represents a window average for a linear depth
profile, and the SSURGO prior side is a genuine thickness-weighted
window mean.

## Usage

``` r
fetch_solus_low_pred_high(aoi_vect, solus_variable, top_depth, bottom_depth)
```

## Arguments

- aoi_vect:

  A
  [`terra::SpatVector`](https://rspatial.github.io/terra/reference/SpatVector-class.html)
  AOI.

- solus_variable:

  A
  [`soilDB::fetchSOLUS()`](http://ncss-tech.github.io/soilDB/reference/fetchSOLUS.md)-recognized
  variable name (e.g. `"claytotal"`, `"dbovendry"`, `"ph1to1h2o"`).

- top_depth, bottom_depth:

  Numeric depth window bounds in cm.

## Value

`list(pred=, low=, high=)`, each a single-layer
[`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
(the window depth-average) or `NULL` if that output type wasn't returned
for every required slice.

## Details

The same window weights are applied to the `prediction` raster and both
`95% ... prediction interval` rasters. Averaging the interval bounds
treats the prediction interval as perfectly depth-correlated
(comonotone) - a mildly conservative (wider) choice for a likelihood,
and the opposite assumption (depth-independence, narrow by ~sqrt(n)) is
equally wrong the other way.
