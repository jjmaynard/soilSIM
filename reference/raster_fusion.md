# Raster-Native Bayesian Fusion (Prior x Likelihood, per Cell)

Per-cell Bayesian fusion of a prior belief distribution (e.g. an
SSURGO-derived percentile raster set) with a likelihood (e.g. an
independent percentile raster set from another source), across a whole
[`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
AOI at once. This is the raster-native counterpart of
`R/bayesian-updating.R` (itself "intentionally NOT called anywhere in
`monte-carlo.R`" - a standalone toolkit), ported from
`code_ref/reanalysis-platform/{bayes_fuse.R, property_fusion_dispatch.R, reanalysis-platform-fusion.R}`.

Reuses `bayesian-updating.R`'s existing scalar/vector fusion functions
directly rather than duplicating them -
[`bayes_update_normal_normal()`](https://jjmaynard.github.io/soilSIM/reference/bayes_update_normal_normal.md),
[`fuse_beta()`](https://jjmaynard.github.io/soilSIM/reference/fuse_beta.md),
[`fuse_gamma()`](https://jjmaynard.github.io/soilSIM/reference/fuse_gamma.md),
[`moments_to_gamma()`](https://jjmaynard.github.io/soilSIM/reference/moments_to_gamma.md)/[`moments_to_beta()`](https://jjmaynard.github.io/soilSIM/reference/moments_to_beta.md)/
[`beta_to_moments()`](https://jjmaynard.github.io/soilSIM/reference/beta_to_moments.md)/[`gamma_to_moments()`](https://jjmaynard.github.io/soilSIM/reference/gamma_to_moments.md),
[`normal_to_lognormal_params()`](https://jjmaynard.github.io/soilSIM/reference/normal_to_lognormal_params.md)/
[`lognormal_to_normal_params()`](https://jjmaynard.github.io/soilSIM/reference/lognormal_to_normal_params.md),
[`bayesian_update()`](https://jjmaynard.github.io/soilSIM/reference/bayesian_update.md),
and `R/distributions.R`'s
[`estimate_ilr_moments_mc()`](https://jjmaynard.github.io/soilSIM/reference/estimate_ilr_moments_mc.md)/[`ilr_inverse()`](https://jjmaynard.github.io/soilSIM/reference/ilr_inverse.md),
plus
[`fuse_bivariate_normal()`](https://jjmaynard.github.io/soilSIM/reference/fuse_bivariate_normal.md)
(also in `bayesian-updating.R`) - are all pure elementwise arithmetic,
so they already work unchanged on `SpatRaster` inputs (confirmed by a
dedicated smoke test before this file was written; see
[`bayes_update_normal_normal()`](https://jjmaynard.github.io/soilSIM/reference/bayes_update_normal_normal.md)'s
own doc comment upstream, which claims exactly this property).

## Deliberately out of scope

[`run_stage1_fusion()`](https://jjmaynard.github.io/soilSIM/reference/run_stage1_fusion.md)/[`run_stage1_fusion_group()`](https://jjmaynard.github.io/soilSIM/reference/run_stage1_fusion_group.md)
(the original top-level orchestrators) were initially **not** ported,
since at that point their dependencies -
[`build_cache_key()`](https://jjmaynard.github.io/soilSIM/reference/build_cache_key.md),
[`cache_get()`](https://jjmaynard.github.io/soilSIM/reference/cache_get.md)/[`cache_set()`](https://jjmaynard.github.io/soilSIM/reference/cache_set.md),
[`fetch_ssurgo_percentiles()`](https://jjmaynard.github.io/soilSIM/reference/fetch_ssurgo_percentiles.md),
[`fetch_solus_percentiles()`](https://jjmaynard.github.io/soilSIM/reference/fetch_solus_percentiles.md)
(itself calling the live
[`soilDB::fetchSOLUS()`](http://ncss-tech.github.io/soilDB/reference/fetchSOLUS.md)
network API) - didn't exist anywhere in this repo. The source bundle's
own `HANDOFF_NOTES.md` confirms it was built as a handoff package for a
*different, sibling project* ("reanalysis-platform"), not soilSIM, with
these exact fetch/cache functions listed as unresolved integration gaps
for that project. Those gaps have since been closed (`R/raster-cache.R`,
`R/ssurgo-simulation.R`, `R/solus-simulation.R`), and
[`run_stage1_fusion()`](https://jjmaynard.github.io/soilSIM/reference/run_stage1_fusion.md)/
[`run_stage1_fusion_group()`](https://jjmaynard.github.io/soilSIM/reference/run_stage1_fusion_group.md)
are now fully implemented below, wiring everything together.
[`fuse_property_adaptive()`](https://jjmaynard.github.io/soilSIM/reference/fuse_property_adaptive.md)
remains the lower-level entry point for callers who already have their
own pre-fetched `prior_value_rasters`/`lik_value_rasters` and want to
skip the fetch-and-cache wrapper.

Also not ported: the global `PROPERTIES` config-list registry from the
source bundle's `config.R`.
[`group_members()`](https://jjmaynard.github.io/soilSIM/reference/group_members.md)/[`fuse_property_adaptive()`](https://jjmaynard.github.io/soilSIM/reference/fuse_property_adaptive.md)
below take a per-call config list/member vector directly instead,
reusing soilSIM's own existing `config$monte_carlo$composition_groups`
convention (see `R/distributions.R`'s
[`resolve_composition_groups()`](https://jjmaynard.github.io/soilSIM/reference/resolve_composition_groups.md))
rather than introducing a second, parallel config schema.
