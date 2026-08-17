# Simulate values from a distribution reconstructed from summary percentiles

Dispatches to one of several methods for turning a handful of known
percentiles into simulated draws. See the file-level `@description`
above for how the methods differ and when to prefer one over another.

## Usage

``` r
simulate_from_percentiles(
  quantile_df,
  method = c("linear_cdf", "spline", "kde", "beta", "normal"),
  percentile_cols = c("P0", "P5", "P50", "P95", "P100"),
  n = 1000,
  bounds = NULL,
  ...
)
```

## Arguments

- quantile_df:

  A data frame with at least one row and percentile-named columns (e.g.
  "P0", "P5", "P50", "P95", "P100").

- method:

  One of "linear_cdf" (default), "spline", "kde", "beta", "normal".

- percentile_cols:

  Character vector of candidate percentile column names.

- n:

  Number of samples to draw.

- bounds:

  Optional length-2 vector of (lower, upper) bounds. Required for method
  = "beta"; used as padding bounds for "spline" if given.

- ...:

  Additional method-specific arguments (`sample_size` for "kde").

## Value

A numeric vector of `n` simulated values.
