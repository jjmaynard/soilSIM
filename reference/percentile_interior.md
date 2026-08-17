# Interior percentiles only (excludes p=0/p=1, where the metalog logit transform is undefined) - matches `fit_metalog_linear_raster()`'s own restriction.

Interior percentiles only (excludes p=0/p=1, where the metalog logit
transform is undefined) - matches
[`fit_metalog_linear_raster()`](https://jjmaynard.github.io/soilSIM/reference/fit_metalog_linear_raster.md)'s
own restriction.

## Usage

``` r
percentile_interior(probs)
```

## Arguments

- probs:

  Numeric probability vector.

## Value

`list(idx=, probs=)`.
