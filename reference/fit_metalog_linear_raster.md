# Fit a metalog distribution via linear solve, vectorized as raster arithmetic.

For a standard percentile setup where the number of interior percentiles
(excluding p=0/p=1, where logit is undefined) equals the number of
metalog terms, the fit is an EXACTLY-DETERMINED linear system
(`solve(Y, z)`), not an optimization - same closed-form math as
`R/distributions.R`'s
[`fit_metalog_linear()`](https://jjmaynard.github.io/soilSIM/reference/fit_metalog_linear.md),
vectorized here as raster arithmetic: `Y`/`Y_inv` are computed once
(fixed, non-spatial probability grid, via the existing
[`metalog_basis_matrix()`](https://jjmaynard.github.io/soilSIM/reference/metalog_basis_matrix.md)),
then applied to different right-hand-side rasters per coefficient.

## Usage

``` r
fit_metalog_linear_raster(value_rasters, probs, bounds, boundedness)

quantile_metalog_linear_raster(fit, q, bounds, boundedness)
```

## Arguments

- value_rasters:

  List of INTERIOR percentile-value rasters (excluding p=0/p=1 if
  present).

- probs:

  Matching interior probabilities.

- bounds:

  `c(lower, upper)`; required unless `boundedness = "u"`.

- boundedness:

  One of `"u"`/`"sl"`/`"su"`/`"b"` - see `metalog_to_z()`.

- fit:

  Output of `fit_metalog_linear_raster()`.

- q:

  Target quantile probability.

## Details

CAVEAT (from upstream validation): this is the same fast path `rmetalog`
itself takes when its solution is already feasible (implied density
non-negative); `rmetalog` also has an LP-based feasibility-correction
fallback for when it isn't - not reproduced here. Use
[`check_metalog_feasibility_raster()`](https://jjmaynard.github.io/soilSIM/reference/check_metalog_feasibility_raster.md) +
[`quantile_metalog_linear_with_fallback()`](https://jjmaynard.github.io/soilSIM/reference/quantile_metalog_linear_with_fallback.md)
below rather than trusting this fit unconditionally.
