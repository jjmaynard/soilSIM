# Run Stage 1 Fusion for a Whole Compositional Group Jointly

Fuses a compositional group's members (currently only `"texture"`:
clay/sand/silt) jointly via
[`fuse_texture_group()`](https://jjmaynard.github.io/soilSIM/reference/fuse_texture_group.md),
instead of one independent
[`run_fusion()`](https://jjmaynard.github.io/soilSIM/reference/run_fusion.md)
call per member - independent fusion measurably breaks sum-to-100. Each
member's SSURGO/SOLUS percentiles are still fetched/cached per-property
(so a later request for just one member's raw percentiles still hits
cache), but the fusion itself runs once for the whole group, under a
`"texture_group"` cache kind (`"texture_group_raw_draws"` when
`prior_fusion_method = "raw_draws"` was requested - see the section
below for why these must be distinct kinds), with each member's
resulting posterior also seeded into a per-property
`"posterior"`/`"posterior_raw_draws"` cache kind - so three sequential
[`run_fusion()`](https://jjmaynard.github.io/soilSIM/reference/run_fusion.md)
calls for clay, then sand, then silt trigger the joint fetch+fusion
exactly once, not three times.

## Usage

``` r
run_fusion_group(
  aoi_vect,
  group,
  composition_groups,
  property_configs,
  top_depth,
  bottom_depth,
  parallel = FALSE,
  n_cores = NULL,
  seed = NULL
)
```

## Arguments

- aoi_vect:

  A
  [`terra::SpatVector`](https://rspatial.github.io/terra/reference/SpatVector-class.html)
  AOI.

- group:

  A composition group name (e.g. `"texture"`).

- composition_groups:

  A `config$monte_carlo$composition_groups`-shaped list (see
  [`group_members()`](https://jjmaynard.github.io/soilSIM/reference/group_members.md))
  giving this group's member ids in order.

- property_configs:

  A named list, keyed by each member id `composition_groups` lists, of
  per-member config lists (each needs at least `id`/`solus_variable`) -
  the per-call replacement for the original bundle's global `PROPERTIES`
  registry.

- top_depth, bottom_depth:

  Numeric depth bounds in cm.

- parallel, n_cores:

  Passed through to
  [`simulate_ssurgo_mapunit_draws()`](https://jjmaynard.github.io/soilSIM/reference/simulate_ssurgo_mapunit_draws.md)'s
  `parallel`/`n_cores` for the shared draws computation below (only
  relevant when at least one member isn't already disk-cached). Default
  `parallel = FALSE` matches prior behavior exactly.

- seed:

  Optional integer for opt-in determinism, forwarded to
  [`simulate_ssurgo_mapunit_draws()`](https://jjmaynard.github.io/soilSIM/reference/simulate_ssurgo_mapunit_draws.md)
  (and [`set.seed()`](https://rdrr.io/r/base/Random.html) once up
  front). `NULL` (default) = current stochastic behavior.

## Value

The full group result: a named list keyed by member id, each element
`list(posterior = list(value=, ilr_mu=, ilr_Sigma=, percentiles=), dist = "texture_ilr", route = "closed_form_ilr_group", route_detail = NULL, n_fallback_cells = 0)`
(see
[`fuse_texture_group()`](https://jjmaynard.github.io/soilSIM/reference/fuse_texture_group.md)'s
docs - NOT the uniform `(mu,sigma)`/`(alpha,beta)` contract
non-compositional properties get) - or `NULL` if any member's
SSURGO/SOLUS fetch failed.
[`run_fusion()`](https://jjmaynard.github.io/soilSIM/reference/run_fusion.md)'s
own dispatch slices this down to the single requested member
(`group_result[[property_config$id]]`) to keep its own per-property
return contract consistent regardless of `dist`. `percentiles` is THIS
member's own fraction's percentiles (a named list of `SpatRaster`s) -
see
[`fuse_texture_group()`](https://jjmaynard.github.io/soilSIM/reference/fuse_texture_group.md)'s
`@section Sum-to-100 caveat`: unlike `value`, these do NOT generally sum
to 100 across the group's three members.

## Higher-fidelity fusion (opt-in)

Set `prior_fusion_method = "raw_draws"` on ANY member's
`property_configs` entry to fuse the whole group against the real joint
(clay, sand, silt) Monte Carlo draws instead of independent
percentile-reconstructed marginals - see
[`fuse_texture_group()`](https://jjmaynard.github.io/soilSIM/reference/fuse_texture_group.md)'s
docs. It is checked across all members (not just one), since the group
is always fused jointly regardless of which member's config carries the
setting. This forces the shared simulation to run even when every
member's own `"ssurgo"` percentile cache is already warm (the draws that
produced those cached percentiles are not kept - see
[`simulate_ssurgo_mapunit_draws()`](https://jjmaynard.github.io/soilSIM/reference/simulate_ssurgo_mapunit_draws.md)'s
docs). Unset on every member (default) uses the
percentile-reconstruction path; unlike
[`run_fusion()`](https://jjmaynard.github.io/soilSIM/reference/run_fusion.md),
the joint texture-group route does not default to raw draws.
