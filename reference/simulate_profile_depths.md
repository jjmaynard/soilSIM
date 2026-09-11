# Simulate soil profile depths, dispatching on input type

A single entry point for horizon depth/thickness simulation, replacing
the pre-0.2.1 `simulate_profile_depths()` / `_by_collection()` /
`_by_collection_parallel()` trio. Each profile in the collection is
simulated via
[`simulate_and_perturb_soil_profiles()`](https://jjmaynard.github.io/soilSIM/reference/simulate_and_perturb_soil_profiles.md),
then combined into one `SoilProfileCollection`.

## Usage

``` r
simulate_profile_depths(x, ...)

# S3 method for class 'SoilProfileCollection'
simulate_profile_depths(x, seed = 123, parallel = FALSE, n_cores = 6, ...)

# Default S3 method
simulate_profile_depths(x, n_simulations = 100, seed = 123, ...)
```

## Arguments

- x:

  Either a `SoilProfileCollection`
  (`simulate_profile_depths.SoilProfileCollection()`) or a mukey
  string/number (`simulate_profile_depths.default()`, which fetches
  SSURGO data for it and builds the collection before delegating to the
  `SoilProfileCollection` method).

- ...:

  Passed to the method.

- seed:

  An integer to set the random seed for reproducibility (default `123`).

- parallel:

  Logical; run the per-profile simulation in parallel via
  `future`/`future.apply` (default `FALSE`, matching the pre-0.2.1
  default - the parallel path was previously a separate function,
  `simulate_profile_depths(parallel = TRUE)`).

- n_cores:

  Integer, workers to use when `parallel = TRUE` (default `6`).

- n_simulations:

  Integer, the number of triangular draws per component used to derive
  `sim_comppct` via
  [`simulate_component_composition()`](https://jjmaynard.github.io/soilSIM/reference/simulate_component_composition.md)
  (default `100`). `default` method only.

## Value

A `SoilProfileCollection` of simulated/perturbed profiles.

## Note on parallel workers

Dispatched via
[`run_parallel_lapply()`](https://jjmaynard.github.io/soilSIM/reference/run_parallel_lapply.md)
(`R/parallel.R`), this package's shared `future`/`future.apply` helper.
[`future::multisession`](https://future.futureverse.org/reference/multisession.html)
workers are fresh R processes; `future`/`globals` auto-detect that the
per-profile closure calls a `soilSIM`-namespaced function and attach the
package in each worker automatically - this only works when `soilSIM` is
actually installed and attached in the calling session, not merely
`devtools::load_all()`'d. The caller's own
[`future::plan()`](https://future.futureverse.org/reference/plan.html)
(if any) is restored afterward regardless of success/failure.
