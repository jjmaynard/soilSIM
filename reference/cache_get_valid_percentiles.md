# Cache-hit lookup for a `list(values=, probs=)` percentile structure, with shape validation

Combines
[`cache_get()`](https://jjmaynard.github.io/soilSIM/reference/cache_get.md) +
[`unwrap_percentile_list()`](https://jjmaynard.github.io/soilSIM/reference/wrap_percentile_list.md) +
[`is_valid_percentile_list()`](https://jjmaynard.github.io/soilSIM/reference/is_valid_percentile_list.md)
into the one safe operation every
[`fetch_ssurgo_percentiles()`](https://jjmaynard.github.io/soilSIM/reference/fetch_ssurgo_percentiles.md)/[`fetch_solus_percentiles()`](https://jjmaynard.github.io/soilSIM/reference/fetch_solus_percentiles.md)
cache-hit call site needs: a malformed or unreadable hit is treated as a
cache miss (returns `NULL`, so the caller re-fetches) rather than being
trusted as-is - see
[`is_valid_percentile_list()`](https://jjmaynard.github.io/soilSIM/reference/is_valid_percentile_list.md)'s
docs for why that distinction matters.

## Usage

``` r
cache_get_valid_percentiles(key, ttl_seconds = CACHE_TTL_SECONDS)
```

## Arguments

- key, ttl_seconds:

  Passed through to
  [`cache_get()`](https://jjmaynard.github.io/soilSIM/reference/cache_get.md).

## Value

The validated, unwrapped `list(values=, probs=)`, or `NULL` on a
miss/stale/malformed entry.
