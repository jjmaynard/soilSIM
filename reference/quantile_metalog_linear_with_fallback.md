# Metalog quantile with automatic fallback to `linear_cdf` for infeasible cells - validated upstream: zero effect on feasible cells, exact `linear_cdf` match on infeasible ones.

Feasibility is computed ONCE per fit (via
[`check_metalog_feasibility_raster()`](https://jjmaynard.github.io/soilSIM/reference/check_metalog_feasibility_raster.md))
and passed in rather than recomputed per call, since it depends only on
the fitted coefficients, not on which quantile `q` is being evaluated.

## Usage

``` r
quantile_metalog_linear_with_fallback(
  fit,
  infeasible_r,
  full_value_rasters,
  full_probs,
  q,
  bounds,
  boundedness
)
```

## Arguments

- fit:

  Output of
  [`fit_metalog_linear_raster()`](https://jjmaynard.github.io/soilSIM/reference/fit_metalog_linear_raster.md)
  (fit on INTERIOR percentiles).

- infeasible_r:

  Output of
  [`check_metalog_feasibility_raster()`](https://jjmaynard.github.io/soilSIM/reference/check_metalog_feasibility_raster.md)
  for this same `fit`.

- full_value_rasters, full_probs:

  The FULL percentile set (including p=0/p=1 if available) -
  `linear_cdf` has no restriction there and benefits from every
  available knot.

- q:

  Target quantile probability.

- bounds, boundedness:

  Same as
  [`fit_metalog_linear_raster()`](https://jjmaynard.github.io/soilSIM/reference/fit_metalog_linear_raster.md)'s
  arguments.

## Value

list(value = blended quantile raster, used_fallback = infeasible_r
verbatim)
