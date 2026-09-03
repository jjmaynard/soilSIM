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
  which bbox-widens for its tabular by-mukey-list query -
  [`soilDB::mukey.wcs()`](http://ncss-tech.github.io/soilDB/reference/mukey.wcs.md)
  bbox-clips internally regardless of input polygon shape, so both
  resolve to the same extent; passing the footprint here avoids widening
  the *cached* grid's extent unnecessarily for callers that only need
  the AOI itself covered).

## Value

A single-layer, categorical (factor)
[`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
named `"mukey"`, or `NULL` if the grid fetch returns nothing.

## Implementation note

Previously this function *also* issued a separate
`soilDB::SDA_spatialQuery(what="mupolygon")` vector fetch and rasterized
it onto `mukey.wcs()`'s own grid purely as a template, discarding that
grid's own cell values -
[`soilDB::mukey.wcs()`](http://ncss-tech.github.io/soilDB/reference/mukey.wcs.md)'s
own documented example
([`?soilDB::mukey.wcs`](http://ncss-tech.github.io/soilDB/reference/mukey.wcs.md))
reads mukey codes directly off its returned raster's values
(`unique(values(res))`), confirming no separate polygon fetch/rasterize
step is needed. Dropped that redundant fetch (see
`MUKEY_DRAWS_FUSION_IMPROVEMENT_PLAN.md`, task P1.1) - this function now
issues exactly one network call per cache miss instead of two. Also
disk-cached (`raster-cache.R`, kind `"mukey_grid"`, depth-agnostic key
via
[`mukey_grid_cache_key()`](https://jjmaynard.github.io/soilSIM/reference/mukey_grid_cache_key.md))
so repeated calls for the same AOI across separate top-level user calls
hit zero network calls.
