# Calculate Confidence Intervals

Calculates confidence intervals for various statistics.

## Usage

``` r
calculate_confidence_intervals(
  data,
  statistic = "mean",
  confidence_level = 0.95,
  method = "normal",
  n_bootstrap = 1000
)
```

## Arguments

- data:

  Input data

- statistic:

  Statistic to calculate ("mean", "median", "proportion")

- confidence_level:

  Confidence level (default = 0.95)

- method:

  Method for calculation ("normal", "bootstrap", "t")

- n_bootstrap:

  Number of bootstrap samples

## Value

Confidence interval
