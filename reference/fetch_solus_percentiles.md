# Fetch SOLUS100 Percentile-Value Rasters for an AOI

The top-level SOLUS "likelihood" entry point for `R/raster-fusion.R`'s
[`fuse_property_adaptive()`](https://jjmaynard.github.io/soilSIM/reference/fuse_property_adaptive.md).

## Usage

``` r
fetch_solus_percentiles(aoi_vect, solus_variable, top_depth, bottom_depth)
```

## Arguments

- aoi_vect:

  A
  [`terra::SpatVector`](https://rspatial.github.io/terra/reference/SpatVector-class.html)
  AOI.

- solus_variable:

  A
  [`soilDB::fetchSOLUS()`](http://ncss-tech.github.io/soilDB/reference/fetchSOLUS.md)-recognized
  variable name.

- top_depth, bottom_depth:

  Numeric depth window bounds in cm.

## Value

`list(values = list(P025 = <low raster>, P50 = <pred raster>, P975 = <high raster>), probs = c(0.025, 0.5, 0.975))`,
or `NULL` if any of low/pred/high is unavailable.
