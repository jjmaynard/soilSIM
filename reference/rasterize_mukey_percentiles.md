# Merge Per-Mukey Percentile Values onto a Mukey Raster

Merge Per-Mukey Percentile Values onto a Mukey Raster

## Usage

``` r
rasterize_mukey_percentiles(mukey_raster, percentile_by_mukey)
```

## Arguments

- mukey_raster:

  A categorical (factor) mukey
  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  from
  [`fetch_ssurgo_mukey_raster()`](https://jjmaynard.github.io/soilSIM/reference/fetch_ssurgo_mukey_raster.md).

- percentile_by_mukey:

  A data frame with a `mukey` column and one column per percentile (e.g.
  `P05`, `P25`, `P50`, `P75`, `P95`).

## Value

A named list of single-layer
[`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)s,
one per percentile column.
