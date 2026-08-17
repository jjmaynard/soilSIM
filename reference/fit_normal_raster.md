# Fit a Normal distribution's mu/sigma rasters from three percentile-value rasters, then evaluate at a fixed target quantile `q` - pure raster arithmetic throughout.

Fit a Normal distribution's mu/sigma rasters from three percentile-value
rasters, then evaluate at a fixed target quantile `q` - pure raster
arithmetic throughout.

## Usage

``` r
fit_normal_raster(p_lo_r, p50_r, p_hi_r, p_lo, p_hi)

quantile_normal_raster(fit, q)
```

## Arguments

- p_lo_r, p50_r, p_hi_r:

  SpatRasters for the low/median/high percentile values.

- p_lo, p_hi:

  The probabilities `p_lo_r`/`p_hi_r` represent (e.g. 0.025/0.975).

- fit:

  Output of `fit_normal_raster()`.

- q:

  Target quantile probability.

## Value

list(mu = SpatRaster, sigma = SpatRaster)
