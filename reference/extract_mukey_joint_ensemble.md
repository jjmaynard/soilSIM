# Per-Mukey Joint Multivariate Profile Ensemble at Several Depth Windows

The multi-property, multi-depth-window analogue of
[`lookup_mukey_draws`](https://jjmaynard.github.io/soilSIM/reference/lookup_mukey_draws.md)
/
[`lookup_mukey_texture_draws`](https://jjmaynard.github.io/soilSIM/reference/lookup_mukey_texture_draws.md):
runs the (expensive) SSURGO Monte Carlo once via
[`simulate_ssurgo_mapunit_draws`](https://jjmaynard.github.io/soilSIM/reference/simulate_ssurgo_mapunit_draws.md)
and returns, per mukey, the **retained bag of joint realizations** -
every simulated property together, row-aligned across depth windows, so
the KSSL cross-property and cross-depth rank structure is preserved.
This is the input to
[`remarginalize_ensemble_to_posterior()`](https://jjmaynard.github.io/soilSIM/reference/remarginalize_ensemble_to_posterior.md)
(`R/core-fusion.R`), which transforms each pixel's marginals to a fused
posterior while keeping this ensemble's empirical copula.

## Usage

``` r
extract_mukey_joint_ensemble(
  aoi_vect,
  depth_windows,
  n_mc = 1000,
  parallel = FALSE,
  n_cores = NULL,
  config = NULL,
  mukey_raster = NULL,
  draws_by_window = NULL,
  properties = SSURGO_SIM_PROPERTY_COLUMNS,
  requested_properties = NULL,
  seed = NULL
)
```

## Arguments

- aoi_vect:

  A
  [`terra::SpatVector`](https://rspatial.github.io/terra/reference/SpatVector-class.html)
  AOI.

- depth_windows:

  Non-empty list of `c(top, bottom)` numeric pairs (cm, `bottom > top`).

- n_mc, parallel, n_cores, config:

  Passed through to
  [`simulate_ssurgo_mapunit_draws`](https://jjmaynard.github.io/soilSIM/reference/simulate_ssurgo_mapunit_draws.md).

- mukey_raster:

  Optional already-fetched mukey raster for this AOI (cache hit
  otherwise).

- draws_by_window:

  Optional pre-computed
  `simulate_ssurgo_mapunit_draws(..., depth_windows =)` result (named
  list of per-window data frames) - supply to reuse an existing
  simulation or to test offline; when given,
  `aoi_vect`/`n_mc`/`parallel`/`n_cores`/`config` are unused.

- properties:

  Property columns to retain (default `SSURGO_SIM_PROPERTY_COLUMNS`);
  the returned `properties` element is the subset actually present for
  this AOI.

- requested_properties:

  Passed through to
  [`simulate_ssurgo_mapunit_draws()`](https://jjmaynard.github.io/soilSIM/reference/simulate_ssurgo_mapunit_draws.md)
  to restrict the simulation itself (as opposed to `properties`, which
  only filters columns after the fact). `NULL` (default) restricts to
  `properties`; ignored when `draws_by_window` is supplied.

- seed:

  Optional integer for opt-in determinism, forwarded to
  [`simulate_ssurgo_mapunit_draws()`](https://jjmaynard.github.io/soilSIM/reference/simulate_ssurgo_mapunit_draws.md);
  ignored when `draws_by_window` is supplied.

## Value

`list(by_mukey = <named list, keyed by mukey code>, properties = <character>, depth_windows = , window_names = , mukey_raster = )`.
Each `by_mukey` element is
`list(windows = <named list of [n_replicate x n_property] matrices, one per window>, replicate_key = data.frame(cokey, simulation_number))`.
`NULL` if the simulation fails or no mukey has a replicate spanning
every window.

## Details

Replicates are aligned across windows on `(cokey, simulation_number)`; a
replicate that has no thickness overlap with some window (so no
aggregate there) is dropped from that mukey entirely, keeping every
window's matrix the same number of rows.

## See also

[`remarginalize_ensemble_to_posterior()`](https://jjmaynard.github.io/soilSIM/reference/remarginalize_ensemble_to_posterior.md),
[`lookup_mukey_draws`](https://jjmaynard.github.io/soilSIM/reference/lookup_mukey_draws.md)
