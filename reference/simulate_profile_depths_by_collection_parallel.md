# Run Soil Profile Depth Simulations in Parallel for a SoilProfileCollection

Simulates soil profile depths for each profile in a
`SoilProfileCollection` in parallel using the `future`/`future.apply`
packages, then combines the results into a single
`SoilProfileCollection`.

## Usage

``` r
simulate_profile_depths_by_collection_parallel(
  soil_collection,
  seed = 123,
  n_cores = 6
)
```

## Arguments

- soil_collection:

  A SoilProfileCollection object containing soil profile data.

- seed:

  An integer value to set the random seed for reproducibility (default
  is 123).

- n_cores:

  Integer, number of cores to use for parallel processing (default is
  6).

## Value

A SoilProfileCollection object containing the simulated and perturbed
soil profiles.

## Note on parallel workers

Dispatched via
[`run_parallel_lapply()`](https://jjmaynard.github.io/soilSIM/reference/run_parallel_lapply.md)
(`R/parallel-utils.R`), this package's shared `future`/`future.apply`
helper.
[`future::multisession`](https://future.futureverse.org/reference/multisession.html)
workers are fresh R processes; `future`/`globals` auto-detect that
`simulate_single_profile()` calls a `soilSIM`-namespaced function and
attach the package in each worker automatically - this only works when
`soilSIM` is actually installed and attached in the calling session, not
merely `devtools::load_all()`'d. The caller's own
[`future::plan()`](https://future.futureverse.org/reference/plan.html)
(if any) is restored afterward regardless of success/failure.

## Examples

``` r
if (FALSE) { # \dontrun{
  simulated_profiles <- simulate_profile_depths_by_collection_parallel(
    soil_collection, seed = 123, n_cores = 6
  )
  print(simulated_profiles)
} # }
```
