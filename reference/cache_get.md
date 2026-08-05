# Retrieve a cached value if present and not older than `ttl_seconds`

Retrieve a cached value if present and not older than `ttl_seconds`

## Usage

``` r
cache_get(key, ttl_seconds = CACHE_TTL_SECONDS)
```

## Arguments

- key:

  A cache key from
  [`build_cache_key()`](https://jjmaynard.github.io/soilSIM/reference/build_cache_key.md).

- ttl_seconds:

  Maximum cache age in seconds before a hit is treated as stale.

## Value

The cached value, or `NULL` on a cache miss/stale entry/read failure.
