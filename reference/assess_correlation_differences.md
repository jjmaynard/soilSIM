# Assess Correlation Differences

Unlike most of its neighbors in this section, this one performs real
computation: element-wise absolute difference between two correlation
matrices, summarized as max/mean/RMSE.

## Usage

``` r
assess_correlation_differences(orig_cor, sim_cor, criteria)
```

## Arguments

- orig_cor:

  Original correlation matrix.

- sim_cor:

  Simulated correlation matrix (same dimensions as `orig_cor`).

- criteria:

  Unused.

## Value

List with `max_difference`, `mean_difference`, `correlation_rmse`.
