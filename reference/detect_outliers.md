# Detect Outliers

Outlier detection using multiple methods.

## Usage

``` r
detect_outliers(data, method = "iqr", threshold = NULL, return_indices = FALSE)
```

## Arguments

- data:

  Input data

- method:

  Detection method ("iqr", "zscore", "modified_zscore")

- threshold:

  Threshold parameter

- return_indices:

  Whether to return indices instead of logical vector

## Value

Outlier indicators or indices
