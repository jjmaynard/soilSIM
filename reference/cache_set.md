# Store a value in the cache under `key`

Store a value in the cache under `key`

## Usage

``` r
cache_set(key, kind, value)
```

## Arguments

- key:

  A cache key from
  [`build_cache_key()`](https://jjmaynard.github.io/soilSIM/reference/build_cache_key.md).

- kind:

  Unused beyond documenting intent at call sites (matches the original
  bundle's 3-argument `cache_set(key, kind, value)` calling
  convention) - the value is stored keyed only by `key`, since
  [`build_cache_key()`](https://jjmaynard.github.io/soilSIM/reference/build_cache_key.md)
  already encodes `kind`.

- value:

  The R object to cache (anything
  [`saveRDS()`](https://rdrr.io/r/base/readRDS.html) can serialize,
  including
  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  objects via their `wrap()`-compatible representation - see Known
  limitation below).

## Value

`TRUE` on success, `FALSE` on a write failure (never errors).

## Known limitation

[`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)/`SpatVector`
objects hold an external pointer to in-memory/on-disk GDAL state that
does not survive a plain
[`saveRDS()`](https://rdrr.io/r/base/readRDS.html)/[`readRDS()`](https://rdrr.io/r/base/readRDS.html)
round-trip as-is. Callers caching raster values must
[`terra::wrap()`](https://rspatial.github.io/terra/reference/wrap.html)
them first and
[`terra::unwrap()`](https://rspatial.github.io/terra/reference/wrap.html)
on read
([`cache_get()`](https://jjmaynard.github.io/soilSIM/reference/cache_get.md)/`cache_set()`
do not do this automatically, since not every cached value is a raster -
e.g. Monte Carlo draw data frames are not).
