# Fuse an Already-Fetched Texture Group (Stage 1 group-fusion tail)

The portion of
[`run_fusion_group()`](https://jjmaynard.github.io/soilSIM/reference/run_fusion_group.md)
that runs once every member's aligned SSURGO prior and SOLUS likelihood
are assembled into `fetched`: optionally prepare the real joint
per-mukey texture draws for raw-draws fusion, then fuse the group via
[`fuse_texture_group()`](https://jjmaynard.github.io/soilSIM/reference/fuse_texture_group.md).
Factored out so
[`run_fusion_group()`](https://jjmaynard.github.io/soilSIM/reference/run_fusion_group.md)
and
[`run_fusion_multiproperty()`](https://jjmaynard.github.io/soilSIM/reference/run_fusion_multiproperty.md)
share one group-fusion code path.

## Usage

``` r
stage1_fuse_texture_group_from_fetched(
  fetched,
  want_raw_draws = FALSE,
  shared_draws = NULL,
  shared_mukey_raster = NULL
)
```

## Arguments

- fetched:

  A list (one entry per group member, in
  [`group_members()`](https://jjmaynard.github.io/soilSIM/reference/group_members.md)
  order) of
  `list(id=, prior=<aligned percentile-value rasters>, prior_probs=, lik=<value rasters>, lik_probs=)`.

- want_raw_draws:

  Whether any member requested `prior_fusion_method = "raw_draws"` - the
  raw-draws texture path only runs when this is `TRUE` (a cold
  percentile cache alone can leave `shared_draws` non-`NULL` without
  raw-draws fusion being requested).

- shared_draws, shared_mukey_raster:

  The one-per-group in-memory draws data frame and its native-grid mukey
  raster (already computed by the caller); both required for the
  raw-draws path.

## Value

[`fuse_texture_group()`](https://jjmaynard.github.io/soilSIM/reference/fuse_texture_group.md)'s
return value - a named list keyed by member id.
