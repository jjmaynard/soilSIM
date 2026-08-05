# Resolve a property's distribution family, refining `dist = "auto"` from the AOI's own percentile skew rather than a fixed per-property label.

Computes a skewness proxy directly from the AOI's own percentile
triplet: `((high-median)-(median-low))/(high-low)`, aggregated once
across the AOI (default `median`) - NOT per-cell, which would create map
discontinuity artifacts at family-boundary cells.

## Usage

``` r
resolve_property_dist(
  property_config,
  prior_value_rasters,
  prior_probs,
  aggregate_fun = stats::median
)
```

## Arguments

- property_config:

  List with `dist` and optionally `bounds`/`auto_skew_threshold`.

- prior_value_rasters, prior_probs:

  Prior percentile-value rasters/probabilities.

- aggregate_fun:

  Function to aggregate the AOI-wide skew proxy (default
  [`stats::median`](https://rdrr.io/r/stats/median.html)).

## Value

`list(dist=, dist_source= one of "config"/"auto", skew_proxy= NA_real_ unless dist_source=="auto")`.

## Details

Deliberately narrow: only ever resolves to "beta" (if `bounds` is
configured), "normal", or "lognormal" - never "metalog" or "gamma",
since distinguishing those from a 3-point proxy isn't reliable. NEVER
overrides an explicit (non-"auto") `dist` value.
