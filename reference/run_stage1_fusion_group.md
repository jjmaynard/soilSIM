# Run Stage 1 Fusion for a Whole Compositional Group Jointly

Fuses a compositional group's members (currently only `"texture"`:
clay/sand/silt) jointly via
[`fuse_texture_group()`](https://jjmaynard.github.io/soilSIM/reference/fuse_texture_group.md),
instead of one independent
[`run_stage1_fusion()`](https://jjmaynard.github.io/soilSIM/reference/run_stage1_fusion.md)
call per member - independent fusion measurably breaks sum-to-100. Each
member's SSURGO/SOLUS percentiles are still fetched/cached per-property
(so a later request for just one member's raw percentiles still hits
cache), but the fusion itself runs once for the whole group, under a
`"texture_group"` cache kind, with each member's resulting posterior
also seeded into a per-property `"posterior"` cache kind - so three
sequential
[`run_stage1_fusion()`](https://jjmaynard.github.io/soilSIM/reference/run_stage1_fusion.md)
calls for clay, then sand, then silt trigger the joint fetch+fusion
exactly once, not three times.

## Usage

``` r
run_stage1_fusion_group(
  aoi_vect,
  group,
  composition_groups,
  property_configs,
  top_depth,
  bottom_depth
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

## Value

The full group result: a named list keyed by member id, each element
`list(posterior = list(value=, ilr_mu=, ilr_Sigma=), dist = "texture_ilr", route = "closed_form_ilr_group", route_detail = NULL, n_fallback_cells = 0)`
(see
[`fuse_texture_group()`](https://jjmaynard.github.io/soilSIM/reference/fuse_texture_group.md)'s
docs - NOT the uniform `(mu,sigma)`/`(alpha,beta)` contract
non-compositional properties get) - or `NULL` if any member's
SSURGO/SOLUS fetch failed.
[`run_stage1_fusion()`](https://jjmaynard.github.io/soilSIM/reference/run_stage1_fusion.md)'s
own dispatch slices this down to the single requested member
(`group_result[[property_config$id]]`) to keep its own per-property
return contract consistent regardless of `dist`.
