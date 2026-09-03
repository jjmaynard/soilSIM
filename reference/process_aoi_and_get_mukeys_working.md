# Process AOI and Get Map Unit Keys (Working Version)

Handles spatial processing exactly as in working code

## Usage

``` r
process_aoi_and_get_mukeys_working(
  aoi_wkt,
  verbose = FALSE,
  mukey_raster = NULL
)
```

## Arguments

- aoi_wkt:

  Well-Known Text representation of area of interest

- verbose:

  Logical; provide progress messages

- mukey_raster:

  Optional, already-fetched
  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  of mukey codes for this same AOI (e.g. from
  [`fetch_ssurgo_mukey_raster`](https://jjmaynard.github.io/soilSIM/reference/fetch_ssurgo_mukey_raster.md)).
  When supplied, this function skips its own
  [`soilDB::mukey.wcs()`](http://ncss-tech.github.io/soilDB/reference/mukey.wcs.md)
  call entirely and derives `mu`/`mukey_list` directly from it - see
  `MUKEY_DRAWS_FUSION_IMPROVEMENT_PLAN.md` task P1.2. `NULL` (default)
  preserves the original behavior exactly.

## Value

List with processed AOI, map units, and mukey list
