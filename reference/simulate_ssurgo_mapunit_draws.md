# Monte Carlo-Simulate SSURGO Property Draws for an AOI

Orchestrates the full SSURGO percentile-prior pipeline for one
AOI/depth-window: fetch (cached) tabular SSURGO data, infill missing
values, derive `genhz`, simulate component composition
([`sim_component_comp()`](https://jjmaynard.github.io/soilSIM/reference/sim_component_comp.md))
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
  n_cores = NULL
)

SSURGO_SIM_PROPERTY_COLUMNS
```

## Format

An object of class `character` of length 10.

## Arguments

- aoi_vect:

  A
  [`terra::SpatVector`](https://rspatial.github.io/terra/reference/SpatVector-class.html)
  AOI.

- top_depth, bottom_depth:

  Numeric depth window bounds in cm.

- n_mc:

  Number of triangular draws
  [`sim_component_comp()`](https://jjmaynard.github.io/soilSIM/reference/sim_component_comp.md)
  uses per component (default 1000).

- parallel, n_cores:

  Passed through to
  [`maybe_adjust_soil_data_depth_trend()`](https://jjmaynard.github.io/soilSIM/reference/maybe_adjust_soil_data_depth_trend.md)'s
  `parallel`/`n_cores` - the depth-trend GP adjustment step is this
  function's dominant cost for AOIs with many cokeys, and each cokey's
  GP fitting is independent of every other cokey's. Default
  `parallel = FALSE` matches prior behavior exactly.

## Value

A data frame, one row per `mukey`/`cokey`/`simulation_number` replicate,
with simulated property columns aggregated over the depth window - or
`NULL` if the tabular fetch fails.
