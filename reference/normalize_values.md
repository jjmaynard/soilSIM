# Normalize Values

Data normalization with multiple methods.

## Usage

``` r
normalize_values(x, method = "minmax", center = TRUE, scale = TRUE)
```

## Arguments

- x:

  Input values

- method:

  Normalization method ("minmax", "zscore", "robust", "quantile")

- center:

  Center values (for zscore and robust)

- scale:

  Scale values (for zscore and robust)

## Value

Normalized values

## Usage note

Fully implemented and exported, but not currently called from anywhere
else in the package (confirmed by source grep) - standalone public API
for callers who need it, not dead/broken code.
