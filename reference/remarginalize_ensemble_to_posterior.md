# Per-Pixel Ensemble Re-Marginalization to a Fused Posterior

Rank-preserving marginal transform: for each property and depth window,
maps every retained source realization through `F_prior` (the mukey
ensemble's empirical CDF) then `F_posterior^-1` (the pixel's fused
posterior, piecewise-linear on its percentile knots). Each realization
keeps its rank position in every property and window, so the ensemble's
empirical copula (cross-property, cross-depth rank correlations) carries
over; only the marginals become the pixel's posterior.

## Usage

``` r
remarginalize_ensemble_to_posterior(
  mukey_ensemble,
  posterior_by_property_window,
  n_out = 250,
  summarize = TRUE,
  probs = c(0.05, 0.25, 0.5, 0.75, 0.95)
)
```

## Arguments

- mukey_ensemble:

  An
  [`extract_mukey_joint_ensemble()`](https://jjmaynard.github.io/soilSIM/reference/extract_mukey_joint_ensemble.md)
  result.

- posterior_by_property_window:

  Named list `[[property]][[window]]` of
  [`run_fusion()`](https://jjmaynard.github.io/soilSIM/reference/run_fusion.md)
  posteriors (each `list(percentiles = , ...)`), on a common grid. Only
  properties present in BOTH this and the ensemble, and windows present
  in both, are produced.

- n_out:

  Max source realizations to carry per pixel (default 250). If a mukey
  has more ensemble replicates, the first `n_out` are used (a fixed
  subsample - same rows across every property/window, so the copula is
  preserved).

- summarize:

  `TRUE` (default) returns per-pixel percentile rasters (`probs`);
  `FALSE` returns the raw `n_kept`-layer realization stacks (for
  downstream use, e.g.
  [`remarginalize_awc()`](https://jjmaynard.github.io/soilSIM/reference/remarginalize_awc.md)).

- probs:

  Percentile probabilities for the `summarize = TRUE` output.

## Value

`list(percentiles = [[property]][[window]] = <named list of P.. SpatRasters>)`
when `summarize`, else
`list(ensemble = [[property]][[window]] = <n_kept-layer SpatRaster>)`.
Plus `properties`, `window_names`, `n_out`, `n_kept`.

## Documented approximations

The copula / correlation structure stays mukey-level (no sub-mukey
spatially-varying dependence); within-depth-window vertical shape comes
from the source realization (window aggregates only - no horizon
push-down here); the copula is assumed invariant under
re-marginalization and is the SSURGO tabular one (KSSL + joint-copula
vertical correlation) - SOLUS carries no joint information, so nothing
is re-estimated; the preserved correlation is **rank (Spearman)**, not
Pearson; composition (comppct) weighting stays mukey-level (already in
the source ensemble's differential replicate counts). Still strictly
better than zonal aggregation. See
`docs/09_multi_source_raster_fusion_pipeline.md` -\> "Statistical
structure of the per-pixel bridge" for how this shapes a
derived-quantity uncertainty band.

## See also

[`extract_mukey_joint_ensemble()`](https://jjmaynard.github.io/soilSIM/reference/extract_mukey_joint_ensemble.md),
[`zonal_distribution_from_posterior()`](https://jjmaynard.github.io/soilSIM/reference/zonal_distribution_from_posterior.md),
[`remarginalize_awc()`](https://jjmaynard.github.io/soilSIM/reference/remarginalize_awc.md)

Other raster-fusion:
[`align_percentile_probs()`](https://jjmaynard.github.io/soilSIM/reference/align_percentile_probs.md),
[`build_cache_key()`](https://jjmaynard.github.io/soilSIM/reference/build_cache_key.md),
[`cache_get()`](https://jjmaynard.github.io/soilSIM/reference/cache_get.md),
[`cache_set()`](https://jjmaynard.github.io/soilSIM/reference/cache_set.md),
[`check_metalog_feasibility_raster()`](https://jjmaynard.github.io/soilSIM/reference/check_metalog_feasibility_raster.md),
[`fit_beta_mle_newton_raster()`](https://jjmaynard.github.io/soilSIM/reference/fit_beta_mle_newton_raster.md),
[`fit_beta_mom_raster()`](https://jjmaynard.github.io/soilSIM/reference/fit_beta_mom_raster.md),
[`fit_gamma_mom_raster()`](https://jjmaynard.github.io/soilSIM/reference/fit_gamma_mom_raster.md),
[`fit_metalog_linear_raster()`](https://jjmaynard.github.io/soilSIM/reference/fit_metalog_linear_raster.md),
[`fit_normal_raster()`](https://jjmaynard.github.io/soilSIM/reference/fit_normal_raster.md),
[`fuse_adaptive()`](https://jjmaynard.github.io/soilSIM/reference/fuse_adaptive.md),
[`percentiles_from_draws()`](https://jjmaynard.github.io/soilSIM/reference/percentiles_from_draws.md),
[`quantile_metalog_linear_with_fallback_raster()`](https://jjmaynard.github.io/soilSIM/reference/quantile_metalog_linear_with_fallback_raster.md),
[`remarginalize_awc()`](https://jjmaynard.github.io/soilSIM/reference/remarginalize_awc.md),
[`resolve_property_dist()`](https://jjmaynard.github.io/soilSIM/reference/resolve_property_dist.md),
[`wrap_nested_rasters()`](https://jjmaynard.github.io/soilSIM/reference/wrap_nested_rasters.md)
