# Run Stage 1 Fusion for One Property/Depth over an AOI

Fetches the (cached) SSURGO prior (`R/ssurgo-simulation.R`'s
[`fetch_ssurgo_percentiles()`](https://jjmaynard.github.io/soilSIM/reference/fetch_ssurgo_percentiles.md))
and SOLUS likelihood (`R/solus-simulation.R`'s
[`fetch_solus_percentiles()`](https://jjmaynard.github.io/soilSIM/reference/fetch_solus_percentiles.md))
percentiles, aligns the SSURGO grid onto the SOLUS grid via
[`terra::resample()`](https://rspatial.github.io/terra/reference/resample.html),
then fuses via
[`fuse_property_adaptive()`](https://jjmaynard.github.io/soilSIM/reference/fuse_property_adaptive.md).

## Usage

``` r
run_stage1_fusion(
  aoi_vect,
  property_config,
  top_depth,
  bottom_depth,
  composition_groups = NULL,
  property_configs = NULL,
  parallel = FALSE,
  n_cores = NULL
)
```

## Arguments

- aoi_vect:

  A
  [`terra::SpatVector`](https://rspatial.github.io/terra/reference/SpatVector-class.html)
  AOI (projected, e.g. EPSG:5070).

- property_config:

  A list with `id` (used as
  [`fetch_ssurgo_percentiles()`](https://jjmaynard.github.io/soilSIM/reference/fetch_ssurgo_percentiles.md)'s
  property id and as the cache key), `solus_variable` (a
  [`soilDB::fetchSOLUS()`](http://ncss-tech.github.io/soilDB/reference/fetchSOLUS.md)-recognized
  variable name), and
  `dist`/`bounds`/`boundedness`/`auto_skew_threshold`/`composition_group`
  as needed by
  [`fuse_property_adaptive()`](https://jjmaynard.github.io/soilSIM/reference/fuse_property_adaptive.md).

- top_depth, bottom_depth:

  Numeric depth bounds in cm.

- composition_groups, property_configs:

  Only needed when `property_config$composition_group` is set - passed
  through to
  [`run_stage1_fusion_group()`](https://jjmaynard.github.io/soilSIM/reference/run_stage1_fusion_group.md)
  (see its docs).

- parallel, n_cores:

  Passed through to
  [`simulate_ssurgo_mapunit_draws()`](https://jjmaynard.github.io/soilSIM/reference/simulate_ssurgo_mapunit_draws.md)'s
  `parallel`/`n_cores` (via
  [`fetch_ssurgo_percentiles()`](https://jjmaynard.github.io/soilSIM/reference/fetch_ssurgo_percentiles.md)
  or
  [`run_stage1_fusion_group()`](https://jjmaynard.github.io/soilSIM/reference/run_stage1_fusion_group.md),
  whichever path this call takes) - the SSURGO adapter's per-cokey
  depth-trend GP fitting is this pipeline's dominant cost for AOIs with
  many cokeys. Default `parallel = FALSE` matches prior behavior
  exactly.

## Value

`list(prior=, likelihood=, posterior=, dist=, dist_source=, skew_proxy=, route=, route_detail=, n_fallback_cells=)`,
or `NULL` if the SSURGO or SOLUS side failed. `posterior`'s shape
depends on `dist` - see
[`fuse_property_adaptive()`](https://jjmaynard.github.io/soilSIM/reference/fuse_property_adaptive.md)'s
docs. As of `MUKEY_DRAWS_FUSION_IMPROVEMENT_PLAN.md` tasks P3.2-P3.5,
`posterior` also carries a `percentiles` field (a named list of
`SpatRaster`s, `"P01"`..`"P99"` by default - see
[FUSE_POSTERIOR_DEFAULT_PROBS](https://jjmaynard.github.io/soilSIM/reference/FUSE_POSTERIOR_DEFAULT_PROBS.md))
regardless of `dist`, for every non-compositional property - this
function doesn't expose its own `posterior_probs` parameter (matching
the existing convention that `n_samples`/`grid_resolution` also aren't
threaded this far down; the underlying `fuse_*()` route's own default
applies).

## Details

Compositional properties (`property_config$composition_group` set) are
detected up front and delegated to
[`run_stage1_fusion_group()`](https://jjmaynard.github.io/soilSIM/reference/run_stage1_fusion_group.md)
instead, since they must be fetched/fit jointly.

## Fusion fidelity - `prior_fusion_method` (default changed at P2.10)

Controls whether fusion runs against each mukey's real per-cell Monte
Carlo draws (`"raw_draws"`) or a percentile-reconstructed approximation
(`"percentile"`) - see
[`fuse_general_kde()`](https://jjmaynard.github.io/soilSIM/reference/fuse_general_kde.md)'s
docs and `MUKEY_DRAWS_FUSION_IMPROVEMENT_PLAN.md` tasks
P2.2/P2.3/P2.6/P2.10. **As of P2.10, leaving
`property_config$prior_fusion_method` unset defaults to `"raw_draws"`**
for `dist %in% c("auto","normal","beta","gamma","lognormal")` - see
[`resolve_want_raw_draws()`](https://jjmaynard.github.io/soilSIM/reference/resolve_want_raw_draws.md)
for the exact rule and the benchmark numbers behind it. Only an explicit
`dist = "metalog"` keeps the pre-P2.10 default -
[`fuse_metalog_adapter()`](https://jjmaynard.github.io/soilSIM/reference/fuse_metalog_adapter.md)
does have a raw-draws fit as of P2.12, but its cost/benefit is
unbenchmarked (mirroring the texture-group route's P2.11 treatment), so
the default stays conservative until that's measured. Set
`property_config$prior_fusion_method = "percentile"` explicitly to opt
back into the original behavior for any `dist` (e.g. to match older
cached/tested output, or avoid the live simulation dependency);
`"raw_draws"` explicitly forces it on even for `dist = "metalog"`, now
genuinely used there (P2.12), not just accepted-and-ignored.

The real draws are held in memory only for the duration of this one
call, never disk-cached (see
[`simulate_ssurgo_mapunit_draws()`](https://jjmaynard.github.io/soilSIM/reference/simulate_ssurgo_mapunit_draws.md)'s
docs for why) - so raw_draws fusion (whether defaulted or explicit)
forces a fresh simulation even when this AOI/property/depth-window's
`"ssurgo"` percentile cache is already warm (the draws that produced
that cached percentile summary weren't kept). That simulation still runs
at most once per call regardless of cache state - it is shared with the
percentile-cache-population step when both are needed in the same call.
