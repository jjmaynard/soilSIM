# Raster-Native Bayesian Fusion (Prior x Likelihood, per Cell)

Per-cell Bayesian fusion of a prior belief distribution (e.g. an
SSURGO-derived percentile raster set) with a likelihood (e.g. an
independent percentile raster set from another source), across a whole
[`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
AOI at once. This is the raster-native counterpart of
`R/bayesian-updating.R`.

Uses `bayesian-updating.R`'s scalar/vector fusion functions directly -
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
[`fuse_bivariate_normal()`](https://jjmaynard.github.io/soilSIM/reference/fuse_bivariate_normal.md) -
are all pure elementwise arithmetic, so they work unchanged on
`SpatRaster` inputs.

## Entry points

[`run_stage1_fusion()`](https://jjmaynard.github.io/soilSIM/reference/run_stage1_fusion.md)/[`run_stage1_fusion_group()`](https://jjmaynard.github.io/soilSIM/reference/run_stage1_fusion_group.md)
are the top-level orchestrators, wiring the fetch-and-cache layer
(`R/raster-cache.R`, `R/ssurgo-simulation.R`, `R/solus-simulation.R`) to
the fusion core.
[`fuse_property_adaptive()`](https://jjmaynard.github.io/soilSIM/reference/fuse_property_adaptive.md)
is the lower-level entry point for callers who already have their own
pre-fetched `prior_value_rasters`/ `lik_value_rasters` and want to skip
the fetch-and-cache wrapper.

## Configuration

[`group_members()`](https://jjmaynard.github.io/soilSIM/reference/group_members.md)/[`fuse_property_adaptive()`](https://jjmaynard.github.io/soilSIM/reference/fuse_property_adaptive.md)
take a per-call config list/member vector directly, reusing soilSIM's
`config$monte_carlo$composition_groups` convention (see
`R/distributions.R`'s
[`resolve_composition_groups()`](https://jjmaynard.github.io/soilSIM/reference/resolve_composition_groups.md)).
