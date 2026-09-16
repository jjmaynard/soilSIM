# Monte Carlo-Simulate SSURGO Property Draws for an AOI

Orchestrates the full SSURGO percentile-prior pipeline for one
AOI/depth-window: fetch (cached) tabular SSURGO data, infill missing
values, derive `genhz`, simulate component composition
([`simulate_component_composition()`](https://jjmaynard.github.io/soilSIM/reference/simulate_component_composition.md))
and join it onto horizons, simulate correlated properties per cokey
([`simulate_cokey_generalized()`](https://jjmaynard.github.io/soilSIM/reference/simulate_cokey_generalized.md),
using the KSSL reference correlation matrices), optionally remove
organic horizons and apply depth-trend GP adjustment, then aggregate to
the requested depth window.

## Usage

``` r
simulate_ssurgo_mapunit_draws(
  aoi_vect,
  top_depth,
  bottom_depth,
  n_mc = 1000,
  parallel = FALSE,
  n_cores = NULL,
  config = NULL,
  mukey_raster = NULL,
  depth_windows = NULL,
  requested_properties = NULL,
  seed = NULL
)

SSURGO_SIM_PROPERTY_COLUMNS
```

## Format

An object of class `character` of length 15.

## Arguments

- aoi_vect:

  A
  [`terra::SpatVector`](https://rspatial.github.io/terra/reference/SpatVector-class.html)
  AOI.

- top_depth, bottom_depth:

  Numeric depth window bounds in cm.

- n_mc:

  Number of triangular draws
  [`simulate_component_composition()`](https://jjmaynard.github.io/soilSIM/reference/simulate_component_composition.md)
  uses per component (default 1000).

- parallel, n_cores:

  Passed through to
  [`maybe_adjust_soil_data_depth_trend()`](https://jjmaynard.github.io/soilSIM/reference/maybe_adjust_soil_data_depth_trend.md)'s
  `parallel`/`n_cores` - the depth-trend GP adjustment step is this
  function's dominant cost for AOIs with many cokeys, and each cokey's
  GP fitting is independent of every other cokey's. Default
  `parallel = FALSE` matches prior behavior exactly.

- config:

  Optional Monte Carlo config, passed through to
  [`maybe_adjust_soil_data_depth_trend()`](https://jjmaynard.github.io/soilSIM/reference/maybe_adjust_soil_data_depth_trend.md)
  -\>
  [`apply_local_gp_adjustments()`](https://jjmaynard.github.io/soilSIM/reference/apply_local_gp_adjustments.md)
  -\>
  [`apply_gp_depth_trends()`](https://jjmaynard.github.io/soilSIM/reference/apply_gp_depth_trends.md).
  `config$monte_carlo$vertical_correlation_method` (default
  `"joint_copula"` default, or `"gp_quantile_retrofit"`) selects this
  top-level entry point's vertical-correlation method; `NULL` (default)
  resolves to `"joint_copula"`, matching
  [`default_monte_carlo_config()`](https://jjmaynard.github.io/soilSIM/reference/default_monte_carlo_config.md)'s
  own default.

- mukey_raster:

  Optional, already-fetched
  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  of mukey codes for this same `aoi_vect` (e.g. from
  [`fetch_ssurgo_mukey_raster`](https://jjmaynard.github.io/soilSIM/reference/fetch_ssurgo_mukey_raster.md)),
  passed through to
  [`download_ssurgo_tabular()`](https://jjmaynard.github.io/soilSIM/reference/download_ssurgo_tabular.md)
  so its tabular fetch reuses it instead of issuing a second,
  independent
  [`soilDB::mukey.wcs()`](http://ncss-tech.github.io/soilDB/reference/mukey.wcs.md)
  call. Only consulted on a `"ssurgo_tabular"` cache miss. `NULL`
  (default): the tabular fetch issues its own call.

- depth_windows:

  Optional list of `c(top, bottom)` numeric pairs. `NULL` (default) -
  one aggregation over `[top_depth, bottom_depth]`, returning a single
  data frame exactly as before (bit-identical). When supplied, the
  (expensive) per-horizon simulation runs once and is then aggregated
  once per window, returning a **named list** of data frames (names
  `"top-bottom"`). `top_depth`/`bottom_depth` are ignored in this mode -
  pass the overall span for them. Built for
  [`extract_mukey_joint_ensemble()`](https://jjmaynard.github.io/soilSIM/reference/extract_mukey_joint_ensemble.md),
  which needs the SOLUS depth windows from one simulation.

- requested_properties:

  `NULL` (default) simulates every recognized property, bit-identical to
  before. Otherwise a character vector (any vocabulary
  [`normalize_requested_properties()`](https://jjmaynard.github.io/soilSIM/reference/normalize_requested_properties.md)
  accepts) restricting
  [`simulate_cokey_generalized()`](https://jjmaynard.github.io/soilSIM/reference/simulate_cokey_generalized.md)
  to just those (plus texture coupling) - the per-cokey GP depth-trend
  fit then runs once per kept property instead of ~8 times, the
  pipeline's dominant saving. A restricted simulation is only
  *law*-equivalent (not bit-identical) to slicing the same columns from
  a full simulation, and only under the default
  `vertical_correlation_method = "joint_copula"`; under
  `"gp_quantile_retrofit"` (not subset-invariant - its rank retrofit
  keys off the first property in the list) `requested_properties` is
  ignored with a warning and all properties are simulated.

- seed:

  Optional integer. `NULL` (default) keeps the current stochastic
  behavior. When set, `set.seed(seed)` runs once at the top (covering
  the whole sequential RNG stream: component composition, per-cokey
  property/texture draws, sequential depth-trend), and the parallel
  depth-trend path uses it as `future.seed`'s L'Ecuyer seed. Determinism
  is **conditional** on identical upstream SSURGO data
  ([`download_ssurgo_tabular()`](https://jjmaynard.github.io/soilSIM/reference/download_ssurgo_tabular.md)
  can change over time) and does **not** make a `requested_properties`
  subset bit-identical to a full simulation.

## Value

With `depth_windows = NULL`: a data frame, one row per
`mukey`/`cokey`/`simulation_number` replicate, with simulated property
columns aggregated over the depth window - or `NULL` if the tabular
fetch fails. With `depth_windows` supplied: a named list of such data
frames, one per window (or `NULL` on fetch failure). Minor components
that
[`download_ssurgo_tabular()`](https://jjmaynard.github.io/soilSIM/reference/download_ssurgo_tabular.md)'s
underlying query would otherwise drop entirely (real `comppct`, zero
`chorizon` rows in SDA) are included transparently here whenever an AOI
sibling let them be recovered - see
[`download_ssurgo_tabular()`](https://jjmaynard.github.io/soilSIM/reference/download_ssurgo_tabular.md)'s
"Component recovery" section.

## Cache invalidation

This function's own disk cache
([`build_cache_key()`](https://jjmaynard.github.io/soilSIM/reference/build_cache_key.md)/[`cache_get()`](https://jjmaynard.github.io/soilSIM/reference/cache_get.md)/[`cache_set()`](https://jjmaynard.github.io/soilSIM/reference/cache_set.md)
via
[`ssurgo_tabular_cache_key()`](https://jjmaynard.github.io/soilSIM/reference/ssurgo_tabular_cache_key.md))
is keyed by `aoi_vect` **only** -
[`download_ssurgo_tabular()`](https://jjmaynard.github.io/soilSIM/reference/download_ssurgo_tabular.md)
takes no depth-window argument and fetches every horizon for every AOI
mukey - and is separate from
[`download_ssurgo_tabular()`](https://jjmaynard.github.io/soilSIM/reference/download_ssurgo_tabular.md)'s
own `cache_dir` parameter (always `NULL` here). One AOI's tabular
download is shared across every depth window/call requested for that
AOI. Clear the cache entry (or the whole cache directory) to force a
fresh fetch.

This function's SIMULATED output is not disk-cached: every call
re-simulates fresh random draws. Raw-draws fusion reuses draws in memory
within one
[`run_fusion()`](https://jjmaynard.github.io/soilSIM/reference/run_fusion.md)
call (computed once, used for both the percentile cache and
[`lookup_mukey_draws()`](https://jjmaynard.github.io/soilSIM/reference/lookup_mukey_draws.md));
see also
[`run_fusion_group()`](https://jjmaynard.github.io/soilSIM/reference/run_fusion_group.md)'s
`shared_draws` pattern.

## See also

Other raster-fusion:
[`align_percentile_probs()`](https://jjmaynard.github.io/soilSIM/reference/align_percentile_probs.md),
[`build_cache_key()`](https://jjmaynard.github.io/soilSIM/reference/build_cache_key.md),
[`cache_get()`](https://jjmaynard.github.io/soilSIM/reference/cache_get.md),
[`cache_set()`](https://jjmaynard.github.io/soilSIM/reference/cache_set.md),
[`check_metalog_feasibility_raster()`](https://jjmaynard.github.io/soilSIM/reference/check_metalog_feasibility_raster.md),
[`closest_solus_depth_slice()`](https://jjmaynard.github.io/soilSIM/reference/closest_solus_depth_slice.md),
[`fetch_solus_low_pred_high()`](https://jjmaynard.github.io/soilSIM/reference/fetch_solus_low_pred_high.md),
[`fetch_solus_percentiles()`](https://jjmaynard.github.io/soilSIM/reference/fetch_solus_percentiles.md),
[`fetch_solus_percentiles_multiproperty()`](https://jjmaynard.github.io/soilSIM/reference/fetch_solus_percentiles_multiproperty.md),
[`fetch_ssurgo_mukey_raster()`](https://jjmaynard.github.io/soilSIM/reference/fetch_ssurgo_mukey_raster.md),
[`fetch_ssurgo_percentiles()`](https://jjmaynard.github.io/soilSIM/reference/fetch_ssurgo_percentiles.md),
[`fit_beta_mle_newton_raster()`](https://jjmaynard.github.io/soilSIM/reference/fit_beta_mle_newton_raster.md),
[`fit_beta_mom_raster()`](https://jjmaynard.github.io/soilSIM/reference/fit_beta_mom_raster.md),
[`fit_gamma_mom_raster()`](https://jjmaynard.github.io/soilSIM/reference/fit_gamma_mom_raster.md),
[`fit_metalog_linear_raster()`](https://jjmaynard.github.io/soilSIM/reference/fit_metalog_linear_raster.md),
[`fit_normal_raster()`](https://jjmaynard.github.io/soilSIM/reference/fit_normal_raster.md),
[`fuse_adaptive()`](https://jjmaynard.github.io/soilSIM/reference/fuse_adaptive.md),
[`percentiles_from_draws()`](https://jjmaynard.github.io/soilSIM/reference/percentiles_from_draws.md),
[`property_to_sim_column()`](https://jjmaynard.github.io/soilSIM/reference/property_to_sim_column.md),
[`quantile_metalog_linear_with_fallback_raster()`](https://jjmaynard.github.io/soilSIM/reference/quantile_metalog_linear_with_fallback_raster.md),
[`remarginalize_awc()`](https://jjmaynard.github.io/soilSIM/reference/remarginalize_awc.md),
[`remarginalize_ensemble_to_posterior()`](https://jjmaynard.github.io/soilSIM/reference/remarginalize_ensemble_to_posterior.md),
[`resolve_property_dist()`](https://jjmaynard.github.io/soilSIM/reference/resolve_property_dist.md),
[`solus_depth_window_weights()`](https://jjmaynard.github.io/soilSIM/reference/solus_depth_window_weights.md),
[`wrap_nested_rasters()`](https://jjmaynard.github.io/soilSIM/reference/wrap_nested_rasters.md)
