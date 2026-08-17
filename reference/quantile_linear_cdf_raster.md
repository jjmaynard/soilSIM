# Evaluate the piecewise-linear inverse-CDF at a fixed quantile `q` across a raster. Because the percentile breakpoints (`probs`) are the same at every cell, only the two value-rasters bracketing `q` and a single scalar interpolation weight are needed - no per-cell branching required since the bracket is resolved once, outside the raster arithmetic.

Evaluate the piecewise-linear inverse-CDF at a fixed quantile `q` across
a raster. Because the percentile breakpoints (`probs`) are the same at
every cell, only the two value-rasters bracketing `q` and a single
scalar interpolation weight are needed - no per-cell branching required
since the bracket is resolved once, outside the raster arithmetic.

## Usage

``` r
quantile_linear_cdf_raster(value_rasters, probs, q)
```

## Arguments

- value_rasters:

  List of SpatRasters, sorted the same way as `probs`.

- probs:

  Numeric probabilities matching `value_rasters`' order.

- q:

  Target quantile probability.
