# Run several percentile-reconstruction methods on the same data and compare them

Draws `n` samples from each requested method and returns both the raw
samples and a summary-statistics table (via
[`calculate_summary_statistics()`](https://jjmaynard.github.io/soilSIM/reference/calculate_summary_statistics.md))
side by side, so the shape/spread of each reconstruction can be compared
by eye or by downstream tests (e.g.
[`stats::ks.test()`](https://rdrr.io/r/stats/ks.test.html) between pairs
of methods).

## Usage

``` r
compare_percentile_methods(
  quantile_df,
  methods = c("linear_cdf", "spline", "kde"),
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

- methods:

  Character vector of methods to run (see `simulate_from_percentiles`).

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

A list with `samples` (named list of numeric vectors, one per method)
and `summary` (data frame of summary statistics, one row per method).
