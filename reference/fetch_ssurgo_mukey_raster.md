# Fetch a Raster of SSURGO Map Unit Keys for an AOI

Fetch a Raster of SSURGO Map Unit Keys for an AOI

## Usage

``` r
fetch_ssurgo_mukey_raster(aoi_vect)
```

## Arguments

- aoi_vect:

  A
  [`terra::SpatVector`](https://rspatial.github.io/terra/reference/SpatVector-class.html)
  (projected, e.g. EPSG:5070) - the AOI footprint itself, not widened to
  its bounding box (unlike `R/ssurgo-acquisition.R`'s
  [`process_aoi_and_get_mukeys_working()`](https://jjmaynard.github.io/soilSIM/reference/process_aoi_and_get_mukeys_working.md),
  which bbox-widens for its tabular by-mukey-list query - inappropriate
  here, where over-fetching/over-rasterizing the AOI's full bbox would
  waste both a `mukey.wcs()` grid and a `SDA_spatialQuery()` polygon
  fetch).

## Value

A single-layer, categorical (factor)
[`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
named `"mukey"`, or `NULL` if either the grid or the polygon fetch
returns nothing (including when
[`soilDB::SDA_spatialQuery()`](http://ncss-tech.github.io/soilDB/reference/SDA_spatialQuery.md)
errors outright, e.g. because its `sf` dependency isn't loadable in this
session - a real, encountered failure mode, not hypothetical).
