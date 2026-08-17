# Method-of-moments Gamma fit from a list of percentile-value rasters

Thin raster wrapper around the existing scalar
[`moments_to_gamma()`](https://jjmaynard.github.io/soilSIM/reference/moments_to_gamma.md)
(`bayesian-updating.R`) - the mean/variance computation is the only
genuinely raster-specific part (combining a *list* of rasters via
[`Reduce()`](https://rdrr.io/r/base/funprog.html)).

## Usage

``` r
fit_gamma_mom_raster(value_rasters)
```

## Arguments

- value_rasters:

  List of percentile-value SpatRasters.

## Value

`list(shape = SpatRaster, rate = SpatRaster)`.
