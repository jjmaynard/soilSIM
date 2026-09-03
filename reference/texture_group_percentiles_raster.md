# Raster-native ILR-posterior percentile sampler for `fuse_texture_group()`

Draws `n_mc` bivariate-Normal ILR-space samples per cell from the fused
`(mu1, mu2, S11, S12, S22)` posterior - the same distribution
`R/distributions.R`'s
[`sample_ilr_posterior()`](https://jjmaynard.github.io/soilSIM/reference/sample_ilr_posterior.md)
draws from in the scalar case, here vectorized raster-wide via
[`terra::app()`](https://rspatial.github.io/terra/reference/app.html) -
applies each cell's own 2x2 Cholesky factor (closed-form:
`u11 = sqrt(S11)`, `u12 = S12/u11`, `u22 = sqrt(S22 - u12^2)`, matching
base R's [`chol()`](https://rdrr.io/r/base/chol.html) upper-triangular
convention exactly),
[`ilr_inverse()`](https://jjmaynard.github.io/soilSIM/reference/ilr_inverse.md)s
the result to `(clay, sand, silt)` triples, then computes
[`quantile()`](https://rdrr.io/r/stats/quantile.html) across draws per
cell per fraction. Unlike
[`bayesian_update()`](https://jjmaynard.github.io/soilSIM/reference/bayesian_update.md)'s
grid-based routes, this posterior has no discretized-grid shortcut -
it's a genuine bivariate Monte Carlo sampler, so `n_mc` is the only
lever on percentile reliability here.

## Usage

``` r
texture_group_percentiles_raster(
  ilr_mu_r,
  ilr_Sigma_r,
  posterior_probs,
  n_mc = 2000,
  max_cells_per_subchunk = 2000
)
```

## Arguments

- ilr_mu_r:

  2-layer SpatRaster (`mu1`, `mu2`) -
  [`fuse_texture_group_batch()`](https://jjmaynard.github.io/soilSIM/reference/fuse_texture_group_batch.md)'s
  `ilr_mu` output.

- ilr_Sigma_r:

  3-layer SpatRaster (`S11`, `S12`, `S22`) - that function's `ilr_Sigma`
  output.

- posterior_probs:

  Numeric probabilities (0-1) at which to report posterior percentiles.

- n_mc:

  Monte Carlo draw count per cell. Default 2000 matches this file's
  existing `texture_n_mc` convention (the separate MC sample size
  [`fuse_texture_group_batch_core()`](https://jjmaynard.github.io/soilSIM/reference/fuse_texture_group_batch.md)
  uses to ESTIMATE the `(mu1,mu2,S11,S12,S22)` moments in the first
  place) - a coincidentally shared default, not the same computation;
  this parameter controls a distinct, later sampling step that draws
  FROM the already-fused posterior distribution.

- max_cells_per_subchunk:

  Same memory rationale as
  [`fuse_texture_group_batch()`](https://jjmaynard.github.io/soilSIM/reference/fuse_texture_group_batch.md)'s
  own parameter - each cell needs `2*n_mc` standard-normal draws plus a
  `3*n_mc`-row ILR-inverse composition matrix, so large
  [`terra::app()`](https://rspatial.github.io/terra/reference/app.html)
  chunks can exhaust memory.

## Value

`list(clay=, sand=, silt=)`, each a named list of percentile
`SpatRaster`s (names `"P1"`, `"P50"`, ... matching `posterior_probs`).

## Sum-to-100 caveat (expected, not a bug)

The point estimate (`posterior$value`, via
[`ilr_inverse()`](https://jjmaynard.github.io/soilSIM/reference/ilr_inverse.md)
on the fused mean) sums to 100 by construction for every cell.
Per-fraction PERCENTILES do **not** generally sum to 100 - percentiles
of correlated marginals don't add linearly (e.g. clay's P95 and sand's
P95 can both be simultaneously "high" draws from different tail regions
of the joint distribution, so their sum can exceed - or their P5s' sum
can fall short of - 100). This is an expected mathematical property of
taking marginal percentiles of a compositional joint distribution, not a
defect.
