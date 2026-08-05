# Fuse a compositional group's members JOINTLY via ILR fusion (`R/distributions.R`'s `estimate_ilr_moments_mc()`/`ilr_inverse()`, `R/bayesian-updating.R`'s `fuse_bivariate_normal()`), rather than independently via `fuse_beta()` per member - independent fusion measurably breaks sum-to-100 (up to 10.5 percentage points on realistic synthetic data, per upstream validation). The raster counterpart of the already-ported scalar `fuse_texture_group_from_triplets()`.

Fuse a compositional group's members JOINTLY via ILR fusion
(`R/distributions.R`'s
[`estimate_ilr_moments_mc()`](https://jjmaynard.github.io/soilSIM/reference/estimate_ilr_moments_mc.md)/[`ilr_inverse()`](https://jjmaynard.github.io/soilSIM/reference/ilr_inverse.md),
`R/bayesian-updating.R`'s
[`fuse_bivariate_normal()`](https://jjmaynard.github.io/soilSIM/reference/fuse_bivariate_normal.md)),
rather than independently via
[`fuse_beta()`](https://jjmaynard.github.io/soilSIM/reference/fuse_beta.md)
per member - independent fusion measurably breaks sum-to-100 (up to 10.5
percentage points on realistic synthetic data, per upstream validation).
The raster counterpart of the already-ported scalar
[`fuse_texture_group_from_triplets()`](https://jjmaynard.github.io/soilSIM/reference/fuse_texture_group_from_triplets.md).

## Usage

``` r
fuse_texture_group(fetched)
```

## Arguments

- fetched:

  A list (one entry per group member, in
  [`group_members()`](https://jjmaynard.github.io/soilSIM/reference/group_members.md)'s
  order) of
  `list(id=, prior=<named list of ALIGNED percentile-value rasters>, prior_probs=, lik=<value rasters>, lik_probs=)`.

## Value

A named list keyed by each member's `id`:
`list(posterior = list(value = <inverse-ILR point-estimate raster for that fraction>, ilr_mu = <2-layer raster, shared across the group>, ilr_Sigma = <3-layer raster (S11,S12,S22), shared>), dist = "texture_ilr", route = "closed_form_ilr_group", route_detail = NULL, n_fallback_cells = 0)`.
To draw posterior samples for a specific fraction/cell, extract that
cell's `ilr_mu`/`ilr_Sigma` and pass to
[`sample_ilr_posterior()`](https://jjmaynard.github.io/soilSIM/reference/sample_ilr_posterior.md).
