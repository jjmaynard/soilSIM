# Fuse a compositional group's members JOINTLY via ILR fusion (`R/core-distributions.R`'s `estimate_ilr_moments_mc()`/`ilr_inverse()`, `R/core-fusion.R`'s `fuse_bivariate_normal()`), rather than independently via `fuse_beta()` per member - independent fusion measurably breaks sum-to-100 (up to 10.5 percentage points on realistic synthetic data). The raster counterpart of the already-ported scalar `fuse_texture_group_from_triplets()`.

Fuse a compositional group's members JOINTLY via ILR fusion
(`R/core-distributions.R`'s
[`estimate_ilr_moments_mc()`](https://jjmaynard.github.io/soilSIM/reference/estimate_ilr_moments_mc.md)/[`ilr_inverse()`](https://jjmaynard.github.io/soilSIM/reference/ilr_inverse.md),
`R/core-fusion.R`'s
[`fuse_bivariate_normal()`](https://jjmaynard.github.io/soilSIM/reference/fuse_bivariate_normal.md)),
rather than independently via
[`fuse_beta()`](https://jjmaynard.github.io/soilSIM/reference/fuse_beta.md)
per member - independent fusion measurably breaks sum-to-100 (up to 10.5
percentage points on realistic synthetic data). The raster counterpart
of the already-ported scalar
[`fuse_texture_group_from_triplets()`](https://jjmaynard.github.io/soilSIM/reference/fuse_texture_group_from_triplets.md).

## Usage

``` r
fuse_texture_group(
  fetched,
  mukey_raster = NULL,
  mukey_texture_draws = NULL,
  posterior_probs = NULL,
  posterior_n_mc = 2000
)
```

## Arguments

- fetched:

  A list (one entry per group member, in
  [`group_members()`](https://jjmaynard.github.io/soilSIM/reference/group_members.md)'s
  order) of
  `list(id=, prior=<named list of ALIGNED percentile-value rasters>, prior_probs=, lik=<value rasters>, lik_probs=)`.

- mukey_raster, mukey_texture_draws:

  Optional - opts into `prior_fusion_method = "raw_draws"` (see
  [`fuse_texture_group_batch_core()`](https://jjmaynard.github.io/soilSIM/reference/fuse_texture_group_batch.md)'s
  docs). `mukey_raster` must already be aligned to the same grid as
  `fetched`'s percentile rasters (nearest-neighbor resampled - it's
  categorical). `mukey_texture_draws` is a
  [`lookup_mukey_texture_draws()`](https://jjmaynard.github.io/soilSIM/reference/lookup_mukey_texture_draws.md)
  result. `NULL` (default) for either uses the percentile-reconstruction
  path.

- posterior_probs:

  Numeric probabilities (0-1) at which to report each fraction's
  posterior percentiles. `NULL` (default) resolves to
  [FUSE_POSTERIOR_DEFAULT_PROBS](https://jjmaynard.github.io/soilSIM/reference/FUSE_POSTERIOR_DEFAULT_PROBS.md).
  Computed via
  [`texture_group_percentiles_raster()`](https://jjmaynard.github.io/soilSIM/reference/texture_group_percentiles_raster.md) -
  see that function's `@section Sum-to-100 caveat` for why, unlike
  `value`, these do NOT generally sum to 100 across the group's three
  members.

- posterior_n_mc:

  Monte Carlo draw count per cell for the percentile sampler (default
  2000). See
  [`texture_group_percentiles_raster()`](https://jjmaynard.github.io/soilSIM/reference/texture_group_percentiles_raster.md)'s
  docs for why this is a separate lever from the moment-estimation MC
  sample size used elsewhere in this function.

## Value

A named list keyed by each member's `id`:
`list(posterior = list(value = <inverse-ILR point-estimate raster for that fraction>, ilr_mu = <2-layer raster, shared across the group>, ilr_Sigma = <3-layer raster (S11,S12,S22), shared>, percentiles = <named list of percentile SpatRasters for THIS member's own fraction - see the sum-to-100 caveat above>), dist = "texture_ilr", route = "closed_form_ilr_group", route_detail = NULL, n_fallback_cells = 0)`.
To draw posterior samples for a specific fraction/cell directly, extract
that cell's `ilr_mu`/`ilr_Sigma` and pass to
[`sample_ilr_posterior()`](https://jjmaynard.github.io/soilSIM/reference/sample_ilr_posterior.md).
