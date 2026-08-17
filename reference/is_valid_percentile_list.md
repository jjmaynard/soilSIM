# Validate a `list(values=, probs=)` percentile structure

[`cache_get()`](https://jjmaynard.github.io/soilSIM/reference/cache_get.md)
returns whatever a cache hit's `.rds` deserializes to without checking
its shape, so a stale or corrupted cache entry (e.g. left behind by an
interrupted run, or written by a since-changed code path) would
otherwise be trusted as-is by every
[`cache_get()`](https://jjmaynard.github.io/soilSIM/reference/cache_get.md)
caller downstream. A malformed hit with an empty/missing `probs` is
especially dangerous here: it doesn't error where it's read, it silently
propagates into
[`align_percentile_probs()`](https://jjmaynard.github.io/soilSIM/reference/align_percentile_probs.md)'s
[`range()`](https://rdrr.io/r/base/range.html) calls (degrading into "no
non-missing arguments to min/max" warnings) and then into
[`fuse_property_adaptive()`](https://jjmaynard.github.io/soilSIM/reference/fuse_property_adaptive.md)'s
`terra::ncell(prior_value_rasters[[1]])` as a cryptic "subscript out of
bounds" error, far from the actual cause. This validator lets cache-hit
call sites treat a malformed hit as a plain cache miss (re-fetch)
instead.

## Usage

``` r
is_valid_percentile_list(x)
```

## Arguments

- x:

  A candidate value read from
  [`cache_get()`](https://jjmaynard.github.io/soilSIM/reference/cache_get.md),
  already
  [`unwrap_percentile_list()`](https://jjmaynard.github.io/soilSIM/reference/wrap_percentile_list.md)-ed.

## Value

`TRUE` if `x` has the expected
`list(values = <non-empty list of SpatRaster>, probs = <numeric vector, same length as values, finite, in [0,1]>)`
shape, else `FALSE`.
