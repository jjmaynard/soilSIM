# Invert a per-pixel posterior CDF at a per-mukey probability, vectorized over the grid

For a scalar-per-mukey probability raster `u_ras`, returns the raster
whose every cell is that cell's posterior value at its `u`, by
piecewise-linear interpolation on the posterior's percentile knots with
`rule = 2` constant extrapolation - the inverse-direction analogue of
[`sim_linear_cdf_batch()`](https://jjmaynard.github.io/soilSIM/reference/sim_linear_cdf_batch.md).

## Usage

``` r
invert_posterior_cdf_raster(pct, post_probs, u_ras)
```

## Arguments

- pct:

  Named list of percentile-value `SpatRaster`s (ascending).

- post_probs:

  Matching ascending probabilities.

- u_ras:

  A `SpatRaster` of probabilities in (0, 1). **May have many layers** -
  each of the single-layer `pct` rasters is recycled across them, so one
  call inverts a whole realization stack in a fixed
  `~4 * (length(post_probs) - 1)` terra ops regardless of `nlyr(u_ras)`.

## Value

A `SpatRaster` of interpolated posterior values, `nlyr(u_ras)` layers.
