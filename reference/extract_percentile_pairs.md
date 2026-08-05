# Parse percentile columns (e.g. "P0","P5","P50") from a one-row data frame

Shared validation/extraction front end used by every method in
[`simulate_from_percentiles()`](https://jjmaynard.github.io/soilSIM/reference/simulate_from_percentiles.md).
Column names must be of the form "P" where is the percentile (0-100).
Columns that are NA or -1 (the sentinel used elsewhere in this project
for "not available") are dropped.

## Usage

``` r
extract_percentile_pairs(quantile_df, percentile_cols)
```

## Arguments

- quantile_df:

  A data frame with at least one row and percentile-named columns.

- percentile_cols:

  Character vector of candidate percentile column names to use.

## Value

A list with `probs` (sorted probabilities, 0-1) and `values` (matching
values).
