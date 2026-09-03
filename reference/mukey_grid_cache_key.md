# Build a cache key for one AOI's SSURGO mukey grid (depth-independent)

A mukey raster/list depends only on the AOI, not on any depth window -
unlike
[`build_cache_key()`](https://jjmaynard.github.io/soilSIM/reference/build_cache_key.md)'s
other callers (percentile-value rasters, tabular horizon data), which
are genuinely depth-window-specific. Rather than changing
[`build_cache_key()`](https://jjmaynard.github.io/soilSIM/reference/build_cache_key.md)'s
signature (every other call site would need a depth argument that means
nothing there), this calls it with a fixed sentinel depth so every
mukey-grid request for the same AOI maps to the same key regardless of
what depth window the caller happens to be working with.

## Usage

``` r
mukey_grid_cache_key(aoi_vect)
```

## Arguments

- aoi_vect:

  A
  [`terra::SpatVector`](https://rspatial.github.io/terra/reference/SpatVector-class.html)
  (single-feature AOI).

## Value

A character string, safe for use as a filename.
