# Fetch SOLUS100 Low/Prediction/High Rasters for One Variable and Depth Window

Fetch SOLUS100 Low/Prediction/High Rasters for One Variable and Depth
Window

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
or `NULL` if that output type wasn't returned.
