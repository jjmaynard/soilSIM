# Fully general fusion path for `fuse_adaptive()`: per cell, draw samples from both sides' percentiles (via `simulate_from_percentiles(method= "linear_cdf")`), fuse via `bayesian_update()`, and moment-match the posterior samples back into the requested family so the output contract matches the closed-form route.

Fully general fusion path for
[`fuse_adaptive()`](https://jjmaynard.github.io/soilSIM/reference/fuse_adaptive.md):
per cell, draw samples from both sides' percentiles (via
`simulate_from_percentiles(method= "linear_cdf")`), fuse via
[`bayesian_update()`](https://jjmaynard.github.io/soilSIM/reference/bayesian_update.md),
and moment-match the posterior samples back into the requested family so
the output contract matches the closed-form route.

## Usage

``` r
fuse_general_kde(
  prior_value_rasters,
  lik_value_rasters,
  percentile_probs,
  family,
  bounds,
  n_samples,
  grid_resolution,
  mukey_raster = NULL,
  mukey_draws = NULL,
  posterior_probs = NULL
)
```

## Arguments

- prior_value_rasters, lik_value_rasters:

  Lists of percentile-value SpatRasters.

- percentile_probs:

  Numeric probabilities matching the value-raster list order.

- family:

  One of "normal", "beta", "gamma".

- bounds:

  Required for `family = "beta"`.

- n_samples:

  Samples drawn per side, per cell.

- grid_resolution:

  Passed to
  [`bayesian_update()`](https://jjmaynard.github.io/soilSIM/reference/bayesian_update.md).
  `NULL` (the caller-facing default from
  [`fuse_adaptive()`](https://jjmaynard.github.io/soilSIM/reference/fuse_adaptive.md))
  resolves to
  [FUSE_GENERAL_KDE_DEFAULT_GRID_RESOLUTION](https://jjmaynard.github.io/soilSIM/reference/FUSE_GENERAL_KDE_DEFAULT_GRID_RESOLUTION.md),
  not
  [`bayesian_update()`](https://jjmaynard.github.io/soilSIM/reference/bayesian_update.md)'s
  own standalone default of `0.01` - see that constant's docs for why.

- mukey_raster, mukey_draws:

  Optional - opts into `prior_fusion_method = "raw_draws"` (see ).
  `mukey_raster` must already be aligned to the same grid as
  `prior_value_rasters`/`lik_value_rasters` (nearest-neighbor resampled,
  since it's categorical - bilinear would fabricate nonsensical mukey
  codes). `mukey_draws` is a
  [`mukey_draws_lookup()`](https://jjmaynard.github.io/soilSIM/reference/mukey_draws_lookup.md)
  result. When both are supplied, the PRIOR side's samples for a cell
  are drawn from that cell's mukey's real Monte Carlo draws (computed
  once per unique mukey, not once per cell) instead of
  [`sim_linear_cdf_batch()`](https://jjmaynard.github.io/soilSIM/reference/sim_linear_cdf_batch.md)'s
  percentile-reconstructed approximation - a genuine fidelity
  improvement, not just a speed one, since the real draws carry
  information the 5-percentile summary lost. The LIKELIHOOD (SOLUS) side
  is unaffected either way, since it has no raw-draws equivalent (see
  that function's own docs). `NULL` (default) for either parameter
  preserves the original percentile-reconstruction behavior for the
  prior side exactly - bit-identical output to before this parameter
  existed.

- posterior_probs:

  Numeric vector of probabilities (0-1) at which to report posterior
  percentiles, computed via
  [`bayesian_update()`](https://jjmaynard.github.io/soilSIM/reference/bayesian_update.md)'s
  grid-based exact percentile machinery (see that function's docs).
  `NULL` (default) resolves to
  [FUSE_POSTERIOR_DEFAULT_PROBS](https://jjmaynard.github.io/soilSIM/reference/FUSE_POSTERIOR_DEFAULT_PROBS.md) -
  unlike
  [`bayesian_update()`](https://jjmaynard.github.io/soilSIM/reference/bayesian_update.md)'s
  own `NULL` default (which means "no percentiles, preserve the original
  plain-vector return"), this is an internal orchestration parameter
  with no external plain-vector consumers, so its `NULL` means "use the
  standard rich default" instead. Also used to compute `mu`/`sigma` (via
  the exact grid-based mean/var, not
  [`mean()`](https://rdrr.io/r/base/mean.html)/[`var()`](https://rdrr.io/r/stats/cor.html)
  on a resampled vector - a strictly more accurate replacement for the
  previous computation, at no extra cost).
