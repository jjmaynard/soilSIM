# Build a cache key for one AOI/property/depth-window/kind combination

Build a cache key for one AOI/property/depth-window/kind combination

## Usage

``` r
build_cache_key(aoi_vect, id, top_depth, bottom_depth, kind)
```

## Arguments

- aoi_vect:

  A
  [`terra::SpatVector`](https://rspatial.github.io/terra/reference/SpatVector-class.html)
  (single-feature AOI).

- id:

  A property id or composition-group name (e.g. `"ph"`, `"texture"`).

- top_depth, bottom_depth:

  Numeric depth bounds in cm.

- kind:

  A short string distinguishing what's cached for this key (e.g.
  `"ssurgo"`, `"solus"`, `"texture_group"`, `"posterior"`).

## Value

A character string, safe for use as a filename.
