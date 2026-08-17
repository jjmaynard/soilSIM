# Cache-hit lookup for a `run_stage1_fusion_group()` result, with shape validation

The nested-group counterpart to
[`cache_get_valid_percentiles()`](https://jjmaynard.github.io/soilSIM/reference/cache_get_valid_percentiles.md) -
see its docs for the general rationale. Combines
[`cache_get()`](https://jjmaynard.github.io/soilSIM/reference/cache_get.md) +
[`unwrap_nested_rasters()`](https://jjmaynard.github.io/soilSIM/reference/wrap_nested_rasters.md) +
[`is_valid_group_result()`](https://jjmaynard.github.io/soilSIM/reference/is_valid_group_result.md).

## Usage

``` r
cache_get_valid_group_result(key, member_ids, ttl_seconds = CACHE_TTL_SECONDS)
```

## Arguments

- key, ttl_seconds:

  Passed through to
  [`cache_get()`](https://jjmaynard.github.io/soilSIM/reference/cache_get.md).

- member_ids:

  Passed through to
  [`is_valid_group_result()`](https://jjmaynard.github.io/soilSIM/reference/is_valid_group_result.md).

## Value

The validated, unwrapped group result, or `NULL` on a
miss/stale/malformed entry.
