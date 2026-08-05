# Wrap/unwrap every `SpatRaster` found anywhere inside an arbitrarily-nested list

A generic counterpart to
[`wrap_percentile_list()`](https://jjmaynard.github.io/soilSIM/reference/wrap_percentile_list.md)/[`unwrap_percentile_list()`](https://jjmaynard.github.io/soilSIM/reference/wrap_percentile_list.md),
for result structures whose shape isn't the fixed
`list(values=, probs=)` percentile form - e.g.
[`run_stage1_fusion()`](https://jjmaynard.github.io/soilSIM/reference/run_stage1_fusion.md)'s
return value, which nests `SpatRaster`s at `prior$values`,
`likelihood$values`, and inside `posterior` (whose own shape varies by
`dist`). Recurses into every list element, replacing each `SpatRaster`
with
[`terra::wrap()`](https://rspatial.github.io/terra/reference/wrap.html)'s
`PackedSpatRaster` representation (or the reverse via
[`terra::unwrap()`](https://rspatial.github.io/terra/reference/wrap.html))
and leaving everything else untouched - directly addresses
[`cache_set()`](https://jjmaynard.github.io/soilSIM/reference/cache_set.md)'s
own `@section Known limitation:` for callers caching a full fusion
result rather than a raw percentile fetch.

## Usage

``` r
wrap_nested_rasters(x)

unwrap_nested_rasters(x)
```

## Arguments

- x:

  A list, arbitrarily nested, that may contain `SpatRaster` objects at
  any depth.

## Value

`x`, with every `SpatRaster` replaced by its wrapped/unwrapped
equivalent.
