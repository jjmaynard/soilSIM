# Compute a metalog fit's mean/sd via numerical quadrature over its own raw quantile function (`mean = integral of Q(p) dp`, `var = integral of Q(p)^2 dp - mean^2`).

Compute a metalog fit's mean/sd via numerical quadrature over its own
raw quantile function (`mean = integral of Q(p) dp`,
`var = integral of Q(p)^2 dp - mean^2`).

## Usage

``` r
metalog_moments_raster(
  fit,
  infeasible_r,
  full_value_rasters,
  full_probs,
  bounds,
  boundedness,
  p_grid = seq(0.001, 0.999, by = 0.005)
)
```

## Arguments

- fit, infeasible_r, bounds, boundedness:

  Same as
  [`check_metalog_feasibility_raster()`](https://jjmaynard.github.io/soilSIM/reference/check_metalog_feasibility_raster.md)'s
  arguments/output for this `fit`.

- full_value_rasters, full_probs:

  The FULL percentile set (including p=0/p=1 if available) - used only
  for the infeasible-cell Normal fallback.

- p_grid:

  Probability grid for the quadrature.

## Known limitation

This is new glue code with no prior validated version (unlike the rest
of this file's math, which was validated upstream against
`fitdistrplus`/ closed-form references) - spot-check against real data
before trusting it in production, per the original source bundle's own
caveat.

Uses
[`quantile_metalog_linear_raster()`](https://jjmaynard.github.io/soilSIM/reference/fit_metalog_linear_raster.md)
(valid for any `p` in `(0,1)`) rather than the `_with_fallback()`
wrapper for the quadrature itself:
[`terra::ifel()`](https://rspatial.github.io/terra/reference/ifelse.html)
isn't lazy, so if the fallback's
[`quantile_linear_cdf_raster()`](https://jjmaynard.github.io/soilSIM/reference/quantile_linear_cdf_raster.md)
were evaluated at every quadrature grid point, it would error whenever
the grid extends outside `full_probs`'s own range. Instead, the
quadrature integrates the raw metalog quantile function over the full
range, and INFEASIBLE cells fall back to
[`fit_normal_raster()`](https://jjmaynard.github.io/soilSIM/reference/fit_normal_raster.md)'s
closed-form moments computed directly from the same percentile triplet.
