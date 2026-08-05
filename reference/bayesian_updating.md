# Bayesian Updating: Conjugate Fusion and General Grid-KDE Fusion

Standalone toolkit for fusing a prior belief distribution with a
likelihood (e.g. an SSURGO-derived prior against field-measured data),
ported from `code_ref/brdf/bayesian_updating.R` and the validated
`code_ref/reanalysis-platform/{bayes_fuse.R, reanalysis-platform-fusion.R, texture_ilr_fusion.R}`
bundle (there raster-native; here plain-vector/ scalar, since the
underlying math is ordinary base-R arithmetic once `terra::` is stripped
out).

**This module is intentionally NOT called anywhere in `monte-carlo.R`.**
It is a standalone, independently-tested set of building blocks for a
future "update an SSURGO-derived prior against observed/field data"
feature - not wired into the main Monte Carlo simulation pipeline in
this pass.

Three tiers, in increasing order of generality/cost:

1.  [`bayes_update_normal_normal()`](https://jjmaynard.github.io/soilSIM/reference/bayes_update_normal_normal.md) -
    exact, closed-form, Normal-only.

2.  [`fuse_beta()`](https://jjmaynard.github.io/soilSIM/reference/fuse_beta.md)/[`fuse_gamma()`](https://jjmaynard.github.io/soilSIM/reference/fuse_gamma.md) -
    exact, closed-form, same-family-only (each family's kernel
    parameters simply ADD - the same principle
    [`bayes_update_normal_normal()`](https://jjmaynard.github.io/soilSIM/reference/bayes_update_normal_normal.md)
    uses for Normal's natural parameters, applied to Beta/Gamma), with a
    `feasible` flag and a documented moment-based fallback route when
    infeasible.

3.  [`bayesian_update()`](https://jjmaynard.github.io/soilSIM/reference/bayesian_update.md) -
    fully general grid-KDE Bayes' rule (any shape, any family, even
    mismatched sides), at the cost of needing raw samples (not
    distributional parameters) on both sides.

[`fuse_bivariate_normal()`](https://jjmaynard.github.io/soilSIM/reference/fuse_bivariate_normal.md)/[`fuse_texture_group_from_triplets()`](https://jjmaynard.github.io/soilSIM/reference/fuse_texture_group_from_triplets.md)
extend tier 1 to the joint 3-part-composition ILR case (clay/sand/silt
texture, in whichever role order `monte-carlo.R`'s
`composition_groups$texture$members` configures - see
`distributions.R`'s ILR section header), since independent per-fraction
fusion (e.g. three separate
[`fuse_beta()`](https://jjmaynard.github.io/soilSIM/reference/fuse_beta.md)
calls) measurably breaks sum-to-100 (documented upstream: up to 10.5
percentage points on realistic synthetic data).
