# Per-unique-mukey closed-form family fit from real Monte Carlo draws, broadcast to a raster

Precomputes a family-native parameter fit (Normal `mu`/`sigma`, Beta
`alpha`/`beta`, or Gamma `shape`/`rate`) once per UNIQUE mukey directly
from that mukey's real simulated draws, then broadcasts the result onto
every cell sharing that mukey via one vectorized categorical lookup
([`terra::subst()`](https://rspatial.github.io/terra/reference/subst.html)) -
not a per-cell computation, unlike
[`fuse_general_kde()`](https://jjmaynard.github.io/soilSIM/reference/fuse_general_kde.md)'s
raw_draws branch (which needs a per-cell
[`density()`](https://rdrr.io/r/stats/density.html) call and so is
meaningfully slower). Verified cheap
(`MUKEY_DRAWS_FUSION_IMPROVEMENT_PLAN.md` tasks P2.6/P2.10): at 100,172
synthetic cells and 150 unique mukeys, the whole precompute+broadcast
pipeline measured ~0.18s total - negligible next to the existing
family-specific fit costs it complements (beta's Newton-Raphson fit
alone measured ~20s at the same cell count), because cost here scales
with mukey CARDINALITY, not cell count.

## Usage

``` r
mukey_draws_closed_form_fit_raster(
  mukey_raster,
  mukey_draws,
  family,
  bounds = NULL
)
```

## Arguments

- mukey_raster:

  Categorical mukey SpatRaster, aligned to the target grid.

- mukey_draws:

  A
  [`mukey_draws_lookup()`](https://jjmaynard.github.io/soilSIM/reference/mukey_draws_lookup.md)
  result - named list keyed by mukey (character), each element a numeric
  vector of real simulated draws for that mukey.

- family:

  One of "normal", "beta", "gamma", "lognormal". `"lognormal"` (added
  `MUKEY_DRAWS_FUSION_IMPROVEMENT_PLAN.md` task P2.12) fits `mu`/`sigma`
  directly from `log(draws)` - the same LOG-SPACE parameterization
  [`fuse_lognormal_adaptive()`](https://jjmaynard.github.io/soilSIM/reference/fuse_lognormal_adaptive.md)'s
  large-AOI branch already uses internally
  ([`normal_to_lognormal_params()`](https://jjmaynard.github.io/soilSIM/reference/normal_to_lognormal_params.md)'s
  output shape), so the merge there is a drop-in replacement of the
  percentile-triplet-derived log-space fit, not a new shape.
  Non-positive draws are dropped before taking
  [`log()`](https://rdrr.io/r/base/Log.html) (a real lognormal draw is
  always positive; any non-positive value reflects upstream noise, not
  signal) - a mukey left with fewer than 2 positive draws falls back to
  `NA` like every other family here.

- bounds:

  Required for `family = "beta"` - draws are scaled to `[0,1]` before
  fitting, to match
  [`fit_beta_mle_newton_raster()`](https://jjmaynard.github.io/soilSIM/reference/fit_beta_mle_newton_raster.md)'s
  own scaled-space convention (see that function's docs) -
  `qfun`/downstream percentile code already expects `alpha`/`beta` in
  that scaled space.

## Value

A named list of `SpatRaster`s (`mu`/`sigma`, `alpha`/`beta`, or
`shape`/`rate`) - `NA` at any cell whose mukey has no `mukey_draws`
entry. No fallback is applied here; the caller merges this against the
percentile-based fit for such cells (matching every other raw_draws
branch's per-cell degrade-to-percentile-reconstruction convention).
