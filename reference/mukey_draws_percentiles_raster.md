# Per-unique-mukey empirical percentiles from real Monte Carlo draws, broadcast to a raster

Companion to
[`mukey_draws_closed_form_fit_raster()`](https://jjmaynard.github.io/soilSIM/reference/mukey_draws_closed_form_fit_raster.md)
for routes that consume percentile VALUES directly rather than a family
parameter fit -
[`fuse_metalog_adapter()`](https://jjmaynard.github.io/soilSIM/reference/fuse_metalog_adapter.md),
whose metalog fit is an exact linear interpolation through fixed
percentile knots ([`solve()`](https://rdrr.io/r/base/solve.html) on a
basis matrix), not a moment/density fit. There is no "metalog fit to raw
draws" analogous to
[`mukey_draws_closed_form_fit_raster()`](https://jjmaynard.github.io/soilSIM/reference/mukey_draws_closed_form_fit_raster.md)'s
normal/beta/gamma/lognormal cases - the real extension for a
percentile-interpolation method is to interpolate through REAL empirical
percentiles
([`stats::quantile()`](https://rdrr.io/r/stats/quantile.html) on the
actual draws) instead of the percentile-reconstruction values. Same
per-unique-mukey precompute +
[`terra::subst()`](https://rspatial.github.io/terra/reference/subst.html)
broadcast pattern as its sibling - cheap, scales with mukey cardinality,
not cell count.

## Usage

``` r
mukey_draws_percentiles_raster(mukey_raster, mukey_draws, probs)
```

## Arguments

- mukey_raster:

  Categorical mukey SpatRaster, aligned to the target grid.

- mukey_draws:

  A
  [`lookup_mukey_draws()`](https://jjmaynard.github.io/soilSIM/reference/lookup_mukey_draws.md)
  result - named list keyed by mukey (character), each element a numeric
  vector of real simulated draws for that mukey.

- probs:

  Numeric probabilities (0-1) at which to compute each mukey's empirical
  quantile - pass the same knot probabilities the caller's
  [`fit_metalog_linear_raster()`](https://jjmaynard.github.io/soilSIM/reference/fit_metalog_linear_raster.md)
  call will use.

## Value

A named list of `SpatRaster`s, one per `probs` entry (named
`paste0("P", round(probs * 100))`, matching this file's other
percentile-list conventions) - `NA` at any cell whose mukey has no
`mukey_draws` entry. No fallback is applied here; the caller merges this
against the percentile-reconstruction values for such cells, same
convention as every other raw_draws branch.
