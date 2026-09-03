# Raw Per-Mukey JOINT Texture Draws, Keyed by Mukey (clay/sand/silt Row-Aligned)

The texture-group counterpart of
[`mukey_draws_lookup`](https://jjmaynard.github.io/soilSIM/reference/mukey_draws_lookup.md):
instead of one property's values, returns the row-aligned
`(clay_total, sand_total, silt_total)` TRIPLES simulated together for
the same cokey/replicate. This preserves the real cross-property
correlation between the three fractions (from
[`simulate_cokey_generalized()`](https://jjmaynard.github.io/soilSIM/reference/simulate_cokey_generalized.md)'s
KSSL texture correlation matrix) - something
[`fuse_texture_group()`](https://jjmaynard.github.io/soilSIM/reference/fuse_texture_group.md)'s
current percentile-reconstruction fusion discards entirely, since it
draws each fraction independently from its own marginal
percentile-derived normal. Used by
[`fuse_texture_group_batch_core()`](https://jjmaynard.github.io/soilSIM/reference/fuse_texture_group_batch.md)'s
`prior_fusion_method = "raw_draws"` opt-in - see
`MUKEY_DRAWS_FUSION_IMPROVEMENT_PLAN.md` task P2.4.

## Usage

``` r
mukey_texture_draws_lookup(draws)
```

## Arguments

- draws:

  A data frame from
  [`simulate_ssurgo_mapunit_draws`](https://jjmaynard.github.io/soilSIM/reference/simulate_ssurgo_mapunit_draws.md),
  for the same AOI/depth window the caller's mukey raster covers.

## Value

A named list keyed by mukey (as character), each element a matrix with
columns `clay_total`/`sand_total`/`silt_total` (one row per
cokey/replicate with a complete texture triplet for that row - rows with
any missing fraction, e.g. a cokey where texture simulation failed for
that replicate, are dropped rather than kept partially). `NULL` if
`draws` doesn't carry all three texture columns at all (e.g. no member
of the AOI's data had texture data).
