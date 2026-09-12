# Run Stage 1 Fusion for Many Properties over an AOI in One Simulation Pass

The multi-property / multi-depth-window counterpart of
[`run_fusion()`](https://jjmaynard.github.io/soilSIM/reference/run_fusion.md):
runs the expensive SSURGO Monte Carlo simulation **once** - covering
every requested property and every depth window via
[`simulate_ssurgo_mapunit_draws()`](https://jjmaynard.github.io/soilSIM/reference/simulate_ssurgo_mapunit_draws.md)'s
`depth_windows` argument - then fuses each property/window against
SOLUS, seeding the same per-property `"ssurgo"`/`"solus"`/`"posterior"`
disk caches
[`run_fusion()`](https://jjmaynard.github.io/soilSIM/reference/run_fusion.md)
reads. A later single-property
[`run_fusion()`](https://jjmaynard.github.io/soilSIM/reference/run_fusion.md)
call for the same AOI/property/window therefore hits cache instead of
re-simulating.

## Usage

``` r
run_fusion_multiproperty(
  aoi_vect,
  property_configs,
  depth_windows,
  composition_groups = NULL,
  n_mc = 1000,
  parallel = FALSE,
  n_cores = NULL,
  seed = NULL,
  verbose = FALSE,
  simplify = FALSE,
  resampling = c("down", "up")
)
```

## Arguments

- aoi_vect:

  A
  [`terra::SpatVector`](https://rspatial.github.io/terra/reference/SpatVector-class.html)
  AOI (projected, e.g. EPSG:5070).

- property_configs:

  A **named** list, keyed by each config's own `id`, of
  `property_config` lists exactly as
  [`run_fusion()`](https://jjmaynard.github.io/soilSIM/reference/run_fusion.md)
  takes them (each needs `id`/`solus_variable`, optionally
  `dist`/`bounds`/`composition_group`/`prior_fusion_method`/...).
  Configs carrying a `composition_group` are fused jointly per
  [`run_fusion_group()`](https://jjmaynard.github.io/soilSIM/reference/run_fusion_group.md);
  **every** member of any referenced group must have its own entry here.

- depth_windows:

  A non-empty list of `c(top, bottom)` numeric pairs (cm,
  `bottom > top`).

- composition_groups:

  A `config$monte_carlo$composition_groups`-shaped list (see
  [`group_members()`](https://jjmaynard.github.io/soilSIM/reference/group_members.md));
  required if any config sets `composition_group`.

- n_mc, parallel, n_cores:

  Passed to
  [`simulate_ssurgo_mapunit_draws()`](https://jjmaynard.github.io/soilSIM/reference/simulate_ssurgo_mapunit_draws.md).

- seed:

  Optional integer for opt-in determinism. `NULL` (default) = current
  stochastic behavior. When set, `set.seed(seed)` runs once up front and
  `seed` is forwarded to the one shared
  [`simulate_ssurgo_mapunit_draws()`](https://jjmaynard.github.io/soilSIM/reference/simulate_ssurgo_mapunit_draws.md)
  call, so the whole batch is reproducible given identical upstream
  SSURGO/SOLUS data.

- verbose:

  Passed to
  [`fuse_property_adaptive()`](https://jjmaynard.github.io/soilSIM/reference/fuse_property_adaptive.md) -
  default `FALSE` (a multi-leaf batch is otherwise noisy with per-route
  messages).

- simplify:

  If `TRUE`, each leaf is reduced to
  `list(percentiles = <named SpatRasters>)` (what per-pixel
  re-marginalization consumes); default `FALSE` returns the full
  [`run_fusion()`](https://jjmaynard.github.io/soilSIM/reference/run_fusion.md)-shaped
  list per leaf.

- resampling:

  Passed through to
  [`align_prior_likelihood()`](https://jjmaynard.github.io/soilSIM/reference/align_prior_likelihood.md)
  for every leaf's prior/likelihood alignment - `"down"` (default)
  matches all prior behavior; `"up"` preserves SSURGO's native grid
  instead. See
  [`align_prior_likelihood()`](https://jjmaynard.github.io/soilSIM/reference/align_prior_likelihood.md)'s
  docs.

## Value

A nested named list `result[[config_id]][["<top>-<bottom>"]]`. Each leaf
is a
[`run_fusion()`](https://jjmaynard.github.io/soilSIM/reference/run_fusion.md)
return list (or, for group members, the `texture_ilr`-shaped list from
[`run_fusion_group()`](https://jjmaynard.github.io/soilSIM/reference/run_fusion_group.md)),
or `NULL` on that leaf's own failure; the whole call returns `NULL` only
if the shared simulation itself fails.

## Details

Calling
[`run_fusion()`](https://jjmaynard.github.io/soilSIM/reference/run_fusion.md)
in a loop over N properties x M windows runs the (dominant-cost)
per-cokey simulation N\*M times, each discarding all but one of the ~9
jointly-simulated properties; this runs it once. Every leaf shares that
one draw set, so - unlike independent
[`run_fusion()`](https://jjmaynard.github.io/soilSIM/reference/run_fusion.md)
calls - their `NA` masks are mutually consistent.

## See also

[`run_fusion()`](https://jjmaynard.github.io/soilSIM/reference/run_fusion.md),
[`run_fusion_group()`](https://jjmaynard.github.io/soilSIM/reference/run_fusion_group.md),
[`extract_mukey_joint_ensemble()`](https://jjmaynard.github.io/soilSIM/reference/extract_mukey_joint_ensemble.md)
