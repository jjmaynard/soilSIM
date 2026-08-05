# Build a piecewise-linear inverse-CDF sampler from percentile columns

Backwards-compatible wrapper around
`simulate_from_percentiles(..., method = "linear_cdf")`.

## Usage

``` r
generate_inverse_cdf_distribution(
  quantile_df,
  percentile_cols = c("P0", "P5", "P50", "P95", "P100"),
  n = 1000
)
```

## Arguments

- quantile_df:

  A data frame with at least one row and percentile-named columns (e.g.
  "P0", "P5", "P50", "P95", "P100").

- percentile_cols:

  Character vector of candidate percentile column names.

- n:

  Number of samples to draw.
