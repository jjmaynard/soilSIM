# Resolve a property's distribution family from its own percentile skew

Computed directly from this row's own l/r/h (no spatial/AOI aggregation
in a tabular context). Deliberately narrow - never resolves to
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
