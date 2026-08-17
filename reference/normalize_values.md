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
