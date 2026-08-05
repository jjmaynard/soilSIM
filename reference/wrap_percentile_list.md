# Wrap/unwrap the `SpatRaster`s inside a `list(values = <list of SpatRaster>, probs = )` percentile structure for caching - see `cache_set()`'s `@section Known limitation:`.

Wrap/unwrap the `SpatRaster`s inside a
`list(values = <list of SpatRaster>, probs = )` percentile structure for
caching - see
[`cache_set()`](https://jjmaynard.github.io/soilSIM/reference/cache_set.md)'s
`@section Known limitation:`.

## Usage

``` r
wrap_percentile_list(x)

unwrap_percentile_list(x)
```

## Arguments

- x:

  A `list(values=, probs=)` structure, as returned by
  [`fetch_ssurgo_percentiles()`](https://jjmaynard.github.io/soilSIM/reference/fetch_ssurgo_percentiles.md)/
  [`fetch_solus_percentiles()`](https://jjmaynard.github.io/soilSIM/reference/fetch_solus_percentiles.md).
