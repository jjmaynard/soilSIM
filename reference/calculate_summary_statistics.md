# Compute summary statistics (mean, SD, CV, percentiles, quartiles) for a numeric vector

Compute summary statistics (mean, SD, CV, percentiles, quartiles) for a
numeric vector

## Usage

``` r
calculate_summary_statistics(data, percentile_probs = seq(0.1, 0.9, by = 0.1))
```

## Arguments

- data:

  A numeric vector.

- percentile_probs:

  Numeric vector of percentile probabilities (0-1) to compute.

## Value

A one-row data frame of summary statistics, including dynamically-named
percentile columns (e.g. P10, P20, ...).
