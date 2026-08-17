# Resolve a property's distribution family from its own percentile skew

Adapted from `code_ref/reanalysis-platform/property_fusion_dispatch.R`'s
[`resolve_property_dist()`](https://jjmaynard.github.io/soilSIM/reference/resolve_property_dist.md):
computed directly from this row's own l/r/h (no spatial/AOI aggregation
needed in a tabular context). Deliberately narrow - never resolves to
`"metalog"` (that stays an explicit config choice), matching the
reference's own documented restriction.

## Usage

``` r
resolve_property_family(l, r, h, bounds = NULL, skew_threshold = 0.15)
```

## Arguments

- l, r, h:

  Low/representative/high values.

- bounds:

  Optional `c(lower, upper)`; when supplied, always resolves to
  `"beta"`.

- skew_threshold:

  Absolute skew-proxy threshold below which `"normal"` is chosen over
  `"lognormal"`.

## Value

`list(family=, skew_proxy=)`.
