# Method-of-moments Beta fit: alpha/beta rasters are pure arithmetic on mean/variance rasters.

Method-of-moments Beta fit: alpha/beta rasters are pure arithmetic on
mean/variance rasters.

## Usage

``` r
fit_beta_mom_raster(mean_r, var_r)

quantile_beta_mom_raster(fit, q)
```

## Arguments

- mean_r, var_r:

  SpatRasters of the rescaled-to-unit-interval mean/variance.

- fit:

  Output of `fit_beta_mom_raster()`.

- q:

  Target quantile probability.
