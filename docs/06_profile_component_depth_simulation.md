# Profile, Component & Depth Simulation

## Overview

This functional group covers two source files that together simulate the
compositional and morphological variability of SSURGO soil components and
profiles. `R/property-simulation.R` simulates component-composition
percentages (`sim_component_comp()`) and per-cokey correlated soil
properties - including compositional sand/silt/clay texture via the
package's ILR core - from low/representative/high (`_l/_r/_h`) triplets
(`simulate_cokey_generalized()`, built on `simulate_correlated_triangular()`).
`R/depth-simulation.R` simulates horizon depths, horizon-thickness
variability, and boundary distinctness for whole `aqp::SoilProfileCollection`
profiles, using OSD (Official Series Description) data to inform boundary
perturbation, with sequential, parallel (`future`/`future.apply`), and
by-mukey orchestration entry points. Both files were ported from
`code_ref/brdf/property_simulation.R` and `code_ref/brdf/depth_simulation.R`
respectively; both are also the two smaller, standalone-horizon-data helpers
`remove_organic_layer()` and `slice_and_aggregate_soil_data()`
(`property-simulation.R`), which are not tied to the compositional-simulation
machinery but shipped from the same legacy file.

## Core Functions

### property-simulation.R

#### `remove_organic_layer(df)`

- **Parameters**: `df` - a data frame of horizon rows with `cokey`, `hzname`,
  `hzdept_r`, `hzdepb_r`, and optionally `hzdept_l`, `hzdepb_l`, `hzdept_h`,
  `hzdepb_h`.
- **Returns**: A data frame with organic horizons removed and remaining
  depths re-anchored to 0 within each `cokey` group.
- **Behavior**: Groups by `cokey`; within each group, drops rows whose
  `hzname` contains a capital "O", shifts the low/high depth columns
  (`hzdep*_l`/`hzdep*_h`) by the first remaining horizon's original
  `hzdept_r` so they stay consistent, zeroes any low/high top-depth value
  that lands exactly on the new top horizon, then recomputes `hzdept_r`/
  `hzdepb_r` (and the low/high pairs, if present and not entirely `NA`) as
  cumulative sums of the original horizon thicknesses - preserving each
  horizon's thickness while removing the organic layer's contribution to
  depth.

#### `slice_and_aggregate_soil_data(df, depth_ranges = list(c(0, 30), c(30, 100)))`

- **Parameters**:
  - `df` - a data frame with `hzdept_r`/`hzdepb_r` columns plus arbitrary
    numeric property columns.
  - `depth_ranges` - a list of length-2 `c(top, bottom)` vectors defining the
    output depth bins.
- **Returns**: A data frame with one row per requested depth range, holding
  the mean of every numeric column across that range (or an all-`NA` row if
  no sliced depth falls in range).
- **Behavior**: Expands every horizon row into individual 1 cm depth
  increments between its `hzdept_r` and `hzdepb_r` (replicating that row's
  property values at each depth), then for each requested `depth_ranges`
  entry filters the exploded data to `Depth >= top & Depth < bottom` and
  takes `colMeans()` (`na.rm = TRUE`). The output row's `hzdept_r` is forced
  to the bin's `top`; `hzdepb_r` is the bin's `bottom`, or one past the
  actual max sliced depth if data ran out before reaching `bottom`. Carries
  `compname` through untransformed if present.

#### `sim_component_comp(data, n_simulations = 1000)`

- **Parameters**:
  - `data` - a data frame with `mukey`, `cokey`, `compname`, `comppct_l`,
    `comppct_r`, `comppct_h` (one or more rows per component; deduplicated
    internally via `dplyr::distinct()`).
  - `n_simulations` - integer, number of triangular draws per component
    (default `1000`).
- **Returns**: A data frame, one row per component (`cokey`), with the
  input columns plus an added `sim_comppct` column.
- **Behavior**: After deduplicating to one row per `cokey` and infilling any
  missing `comppct_l`/`comppct_h` as `comppct_r` -/+ 2, draws `n_simulations`
  triangular samples of `comppct` per component via `tri_dist()` (low =
  `comppct_l`, high = `comppct_h`, mode = `comppct_r`), then sets
  `sim_comppct <- round(sum(draws) / 100)`. This makes `sim_comppct`
  approximately `round(n_simulations * comppct_r / 100)` in expectation
  (e.g. `comppct_r = 30`, `n_simulations = 1000` gives `sim_comppct` ~= 300)
  - a component-percent-weighted realization count, not an independently
    meaningful statistic, ported verbatim from the legacy source. This
    output is per-**component**, but the depth-simulation functions below
    need `sim_comppct` on every **horizon** row - see Known Limitations.

#### `simulate_correlated_triangular(n, params, correlation_matrix, random_seed = NULL)`

- **Parameters**:
  - `n` - integer, number of samples to generate.
  - `params` - a list of `c(a, b, c)` triples, one per distribution, in
    **`(lower, mode, upper)`** order - note this differs from `tri_dist()`'s
    own `(a = lower, b = upper, c = mode)` argument convention; preserved
    exactly as the legacy source defines it.
  - `correlation_matrix` - a square, positive-semi-definite correlation
    matrix, `length(params)` x `length(params)`.
  - `random_seed` - optional integer seed for reproducibility.
- **Returns**: An `n`-row x `length(params)`-column matrix of correlated
  samples.
- **Behavior**: Draws `n * k` independent standard normal values (this is
  mathematically equivalent to `MASS::mvrnorm(n, mu = rep(0, k), Sigma =
  diag(k))`, so no `MASS` dependency is used), correlates them via
  `chol(correlation_matrix)`, converts each correlated-normal column to a
  uniform via `stats::pnorm()`, then inverts the standard triangular-
  distribution CDF analytically (branching on whether `u` falls below or
  above the mode's CDF value) to map each column onto its own `(a, b, c)`
  triangular distribution. When `c == a` for a given parameter, that column
  is filled with the constant `a` (degenerate triangular).

#### `calculate_mode(x)`

- **Parameters**: `x` - a numeric vector.
- **Returns**: The most frequently-occurring value in `x` (via
  `table()`/`which.max()`).
- **Behavior**: Simple statistical-mode helper; used internally by
  `simulate_cokey_generalized()` to collapse a simulated ILR-coordinate
  vector back down to a single representative value for re-use as a
  triangular-distribution mode.

#### `simulate_cokey_generalized(sim_cokey, correlation_matrices, txt_correlation_matrices = NULL)`

- **Parameters**:
  - `sim_cokey` - a data frame of horizon rows for one cokey, with a
    `genhz` column and a `sim_comppct` column (the number of realizations to
    simulate for that row - see `sim_component_comp()`).
  - `correlation_matrices` - a list of correlation matrices keyed by
    `genhz`, with row/column names drawn from (a subset of) `c("db",
    "wr_3b", "wr_15b", "ilr1", "ilr2", "rfv", "ph", "cec", "soc")`.
  - `txt_correlation_matrices` - an optional list of texture correlation
    matrices keyed by `genhz`, in sand/silt/clay row/column order.
- **Returns**: A data frame stacking simulated property values across all
  input rows and all realizations, with `compname`, `mukey`, `cokey`,
  `hzdept_r`, `hzdepb_r`, `simulation_number`, and `unique_id`
  (`"<cokey>-<simulation_number>"`, zero-padded to 2 digits) columns
  attached.
- **Behavior**: For each row of `sim_cokey`, looks up that row's `genhz`-
  keyed local correlation matrix, then builds a parameter list from whichever
  of bulk density (`dbovendry`), 1/3-bar and 15-bar water retention
  (`wthirdbar`, `wfifteenbar`), rock fragment volume (`rfv`), pH
  (`ph1to1h2o`), CEC (`cec7`), and organic matter (`om`) have complete
  `_l/_r/_h` triplets present on that row. If sand/silt/clay triplets are
  all present and `txt_correlation_matrices` was supplied, texture is
  simulated *first* and separately: `simulate_correlated_triangular()` draws
  `sim_comppct` correlated sand/silt/clay values using the texture
  correlation matrix, these are converted to ILR coordinates via
  `ilr_forward()`, and the resulting `ilr1`/`ilr2` vectors are collapsed to
  `(min, calculate_mode(), max)` triples to serve as ordinary
  triangular-distribution parameters for the *main* simulation. All
  recognized parameters (property triplets plus, if present, `ilr1`/`ilr2`)
  are then simulated jointly in a single `simulate_correlated_triangular()`
  call against the row's local (genhz-keyed) correlation matrix subset to
  just the present parameters. Water retention outputs are rescaled `/ 100`;
  if ILR columns were simulated, they are inverted back to `sand_total`/
  `silt_total`/`clay_total` via `ilr_inverse(total = 100)` and the ILR
  columns dropped. Rows with no recognized property triplets are skipped
  (with a `message()`); per-row simulation errors are caught and logged via
  `message()` rather than aborting the whole call. Results across rows are
  `rbind()`-ed into the final data frame.

### depth-simulation.R

#### `get_aws_data_by_mukey(mukeys)`

- **Parameters**: `mukeys` - a vector of string or numeric map unit keys.
- **Returns**: A data frame of SSURGO horizon data for the given mukeys,
  including aggregated rock fragment volume.
- **Behavior**: Builds and submits a SQL query (via `soilDB::SDA_query()`)
  joining `mapunit`/`component`/`chorizon` (inner joins) and `chfrags` (left
  join), pulling component percentages, horizon morphology/depths, texture,
  bulk density, water retention, and rock-fragment-volume fields. Because a
  horizon can have multiple `chfrags` rows (one per fragment-size class),
  results are grouped by `chkey` and `rfv_l`/`rfv_r`/`rfv_h` are summed
  (`na.rm = TRUE`) before being re-merged onto the main result and
  deduplicated. Deliberately does *not* reuse
  `execute_ssurgo_query_working()` (`R/ssurgo-acquisition.R`): that function
  only selects `chf.fragsize_r`, never `chf.fragvol_l/r/h`, so its rock
  fragment path produces no data - delegating to it here would silently
  drop rock fragment volume.

#### `query_osd_distinctness(horizon_data)`

- **Parameters**: `horizon_data` - a data frame with at least `compname` and
  `hzname` columns.
- **Returns**: A data frame with columns `id` (component name), `hzname`,
  `distinctness` (code, e.g. A/C/G/D), `genhz`, and `bound_sd` (the offset
  value derived from the distinctness code).
- **Behavior**: Derives `genhz` for the input horizon data (via
  `aqp::generalizeHz()`) if not already present, renames `compname` to `id`,
  then calls `soilDB::fetchOSD()` for the unique component names to retrieve
  Official Series Description horizon data. Errors if the fetched OSD data
  lacks `distinctness`/`top`/`bottom`/`hzname`/`id` fields. Computes `genhz`
  for the OSD horizons using a distinctness-appropriate numeric-prefix-aware
  regex pattern, identifies any `genhz` values present in the input data but
  missing from the OSD data and appends placeholder rows (`distinctness =
  NA`) for them, calls `infill_missing_distinctness()` to fill remaining
  `NA` distinctness values, then converts distinctness codes to numeric
  offsets via `aqp::hzDistinctnessCodeToOffset()` into `bound_sd`.

#### `infill_missing_distinctness(horizon_data)`

- **Parameters**: `horizon_data` - a data frame with `hzname` and
  `distinctness` columns.
- **Returns**: The same data frame with missing `distinctness` values filled
  in.
- **Behavior**: Derives `genhz` (via `aqp::generalizeHz()`, generalized
  classes `O/A/E/B/C/Cr/R`) if not already present. For each row with a
  missing `distinctness`, looks up a default by `genhz` first, falling back
  to `hzname` if `genhz` isn't in the lookup table; leaves `NA` if neither
  matches. Default distinctness-by-horizon-group table: O = diffuse, A =
  clear, E = clear, B = gradual, C = gradual, Cr = gradual, R = abrupt.

#### `infill_missing_depth_variability(horizon_data)`

- **Parameters**: `horizon_data` - a data frame with `hzdept_r`, `hzdept_l`,
  `hzdept_h`, `hzdepb_r`, `hzdepb_l`, `hzdepb_h` columns.
- **Returns**: The same data frame with missing low/high depth-bound values
  filled in.
- **Behavior**: Row-wise, replaces any missing `hzdept_l`/`hzdepb_l` with
  `hzdept_r - 2`/`hzdepb_r - 2` and any missing `hzdept_h`/`hzdepb_h` with
  `hzdept_r + 2`/`hzdepb_r + 2`, clamped to a minimum of 0 via `pmax()`.

#### `simulate_soil_profile_top_down(horizon_data)`

- **Parameters**: `horizon_data` - a data frame with `hzname`, `hzdepb_l`,
  `hzdepb_r`, `hzdepb_h` columns, one row per horizon in profile order.
- **Returns**: A data frame with `hzname`, `top`, `bottom` - one simulated
  top/bottom depth pair per horizon.
- **Behavior**: Starts the first horizon's `top` at 0. Walks horizons
  top-to-bottom; each subsequent horizon's `top` is set to the previous
  horizon's simulated `bottom`. Each horizon's `bottom` is drawn from
  `tri_dist(n = 1, a = bottom_low, b = bottom_high, c = bottom_mode)`, where
  `bottom_low`/`bottom_high`/`bottom_mode` are the horizon's `hzdepb_l`/
  `hzdepb_h`/`hzdepb_r` clamped so they never fall below the horizon's own
  simulated `top` (preventing an invalid triangular distribution or a
  negative-thickness horizon); the draw is rounded, and if it still comes
  out below the current `top` it is clamped up to `top`.

#### `simulate_soil_profile_bottom_up(horizon_data)`

- **Parameters**: `horizon_data` - a data frame with `hzname`, `hzdept_l`,
  `hzdept_r`, `hzdept_h`, `hzdepb_l`, `hzdepb_r`, `hzdepb_h` columns, one row
  per horizon in profile order.
- **Returns**: A data frame with `hzname`, `top`, `bottom` - one simulated
  top/bottom depth pair per horizon.
- **Behavior**: Walks horizons bottom-to-top. For the deepest horizon, both
  `bottom` (from its own `hzdepb_l/r/h` triangular distribution) and `top`
  (from its own `hzdept_l/r/h`, clamped so `top <= bottom`) are simulated.
  For every horizon above, `bottom` is fixed to the horizon below's already-
  simulated `top` (guaranteeing no gaps), and `top` is drawn from that
  horizon's own `hzdept_l/r/h` triangular distribution, clamped to not
  exceed the fixed `bottom`. After the loop, the uppermost horizon's `top`
  is forced to 0 (adjusting its `bottom` up if that would invert it), and a
  final pass re-checks/enforces that each horizon's `bottom` exactly equals
  the next horizon's `top`.

#### `simulate_soil_profile_thickness(horizon_data, n_simulations = 500)`

- **Parameters**:
  - `horizon_data` - a data frame with `hzname`, `hzdept_r`, `hzdepb_r`,
    `hzdept_l`, `hzdept_h`, `hzdepb_l`, `hzdepb_h` columns.
  - `n_simulations` - integer, number of repeated simulations to run
    (default `500`).
- **Returns**: A data frame with `hzname`, `top`, `bottom` (representative
  depths from the input data) and `thickness_sd` (the standard deviation of
  simulated horizon thickness across all repeated simulations).
- **Behavior**: First calls `infill_missing_depth_variability()` on the
  input. Then, `n_simulations` times, flips a coin (`stats::runif(1) < 0.5`)
  to choose between `simulate_soil_profile_top_down()` and
  `simulate_soil_profile_bottom_up()` for that iteration, records each
  horizon's simulated `Thickness = bottom - top` plus the method and
  simulation index, and stacks all iterations together. Groups the stacked
  results by `hzname` and computes `sd(Thickness)` (rounded to 2 decimals)
  per horizon, then merges that back onto the horizon's representative
  `hzdept_r`/`hzdepb_r` (renamed to `top`/`bottom`) from the original input.
  This is a variability-estimation pass, distinct from the perturbation
  applied downstream by `simulate_and_perturb_soil_profiles()`.

#### `simulate_and_perturb_soil_profiles(soil_profile)`

- **Parameters**: `soil_profile` - a single-profile `aqp::SoilProfileCollection`.
- **Returns**: A `SoilProfileCollection` with `n_simulations` perturbed
  realizations of the input profile (`n_simulations` derived from the
  unique `sim_comppct` value on the profile's horizons - see Known
  Limitations).
- **Behavior**: Subsets the profile's horizons to a fixed working-column set
  (requires `sim_comppct` to already be present) and derives `genhz`. If the
  profile has only one horizon (e.g. an R-only profile), skips perturbation
  entirely and instead returns `n_simulations` verbatim copies of the
  profile, each renamed with a `_sim_<i>` profile-ID suffix. Otherwise:
  (1) calls `simulate_soil_profile_thickness()` to get a per-horizon
  `thickness_sd`; (2) calls `query_osd_distinctness()` and averages
  `bound_sd` by `id`/`genhz` into a lookup table; (3) merges thickness and
  distinctness data together and back onto the profile's horizons; (4) runs
  `aqp::perturb(n = n_simulations, thickness.attr = "thickness_sd")` to
  generate `n_simulations` thickness-perturbed copies of the whole profile
  at once; (5) runs a *second*, per-copy `aqp::perturb(n = 1,
  boundary.attr = "bound_sd")` pass over each of those copies individually
  to additionally perturb horizon boundary depths using the OSD-derived
  distinctness offsets, then recombines them via `aqp::combine()`; (6) as a
  final guard, clamps each simulated horizon's `top`/`bottom` back within
  the original horizon's `hzdept_l`/`hzdepb_h` bounds (via an internal
  `adjust_out_of_range_profiles()` helper) so perturbation cannot push a
  horizon boundary outside its originally observed low/high range.

#### `simulate_profile_depths_by_collection(soil_collection, seed = 123)`

- **Parameters**:
  - `soil_collection` - a multi-profile `aqp::SoilProfileCollection`.
  - `seed` - integer random seed (default `123`).
- **Returns**: A single combined `SoilProfileCollection` of simulated
  profiles for every input profile.
- **Behavior**: Validates the input class, sets the seed, then loops
  sequentially over each profile in the collection, calling
  `simulate_and_perturb_soil_profiles()` on each one-profile subset and
  combining all results via `aqp::combine()`.

#### `simulate_profile_depths_by_collection_parallel(soil_collection, seed = 123, n_cores = 6)`

- **Parameters**:
  - `soil_collection` - a multi-profile `aqp::SoilProfileCollection`.
  - `seed` - integer random seed (default `123`).
  - `n_cores` - integer, number of parallel workers (default `6`).
- **Returns**: A single combined `SoilProfileCollection`, or `NULL` if an
  error occurs during parallel execution.
- **Behavior**: Sets up a `future::multisession` plan with `n_cores`
  workers, validates the input class, sets the seed, then dispatches one
  `simulate_and_perturb_soil_profiles()` call per profile via
  `future.apply::future_lapply()`. Reverts to `future::sequential` after the
  parallel step, combines results via `aqp::combine()`, and wraps the whole
  parallel block in a `tryCatch()` that logs and returns `NULL` on error
  rather than propagating. See the function's "Note on parallel workers"
  section: `future`/`globals` auto-detection of the `soilSIM`-namespaced
  call inside the worker closure only works if `soilSIM` is installed and
  attached in the calling session (not merely `devtools::load_all()`'d).

#### `simulate_profile_depths_by_mukey(mukey, n_simulations = 100, seed = 123)`

- **Parameters**:
  - `mukey` - string or numeric map unit key to query.
  - `n_simulations` - integer (default `100`); **currently unused** - see
    Known Limitations.
  - `seed` - integer random seed (default `123`).
- **Returns**: A combined `SoilProfileCollection` of simulated/perturbed
  profiles for every component in the mukey.
- **Behavior**: Queries `get_aws_data_by_mukey()` for the mukey, errors if
  no data is returned, sets the seed, derives an `id` column from
  `compname` (required because `query_osd_distinctness()` and the internal
  `adjust_out_of_range_profiles()`/`evaluate_simulated_depths()` helpers all
  hardcode a literal `id` column), promotes the data frame to a
  `SoilProfileCollection` via `aqp::depths(mu_data) <- id ~ hzdept_r +
  hzdepb_r` and sets `aqp::hzdesgnname(mu_data) <- "hzname"`, then loops
  over each component profile calling `simulate_and_perturb_soil_profiles()`
  and combining the results via `aqp::combine()`.

#### `evaluate_simulated_depths(simulated_profiles, horizon_data)`

- **Parameters**:
  - `simulated_profiles` - a `SoilProfileCollection` of simulated profiles.
  - `horizon_data` - the original horizon data frame, with `hzname`,
    `hzdept_l`, `hzdepb_h` columns.
- **Returns**: A data frame of just the horizon rows (mukey/cokey/compname/
  hzname/top/bottom) whose simulated depths fall outside the original
  low/high range.
- **Behavior**: Extracts `mukey`/`cokey`/`compname`/`hzname`/`top`
  (`hzdept_r`)/`bottom` (`hzdepb_r`) from the simulated profiles, merges by
  `hzname` with the original data's `hzdept_l`/`hzdepb_h` bounds, flags
  `out_of_range <- (top < hzdept_l) | (bottom > hzdepb_h)`, and returns only
  the flagged rows - a quality-control check on
  `simulate_and_perturb_soil_profiles()`'s output (note that function
  already clamps to this same range internally via
  `adjust_out_of_range_profiles()`, so this is primarily useful for
  auditing/validating that clamp rather than expecting violations in normal
  use).

## Internal Connections

```
property-simulation.R
======================
tri_dist() [R/distributions.R]
    │
    ├──► sim_component_comp() ──► one row per cokey, + sim_comppct
    │
    └──► simulate_correlated_triangular() ◄── ilr_forward()/ilr_inverse() [R/distributions.R]
                  │                                (texture only)
                  ▼
         simulate_cokey_generalized() ──► one row per horizon x realization,
                                           + compname/mukey/cokey/hzdept_r/
                                             hzdepb_r/simulation_number/unique_id

remove_organic_layer() ─────────────────► standalone horizon-cleanup helper
slice_and_aggregate_soil_data() ────────► standalone depth-binning helper
calculate_mode() ────────────────────────► used inside simulate_cokey_generalized()
                                            (collapses simulated ILR draws to a
                                            triangular-mode parameter)


    sim_component_comp() output (per-cokey sim_comppct)
                  │
                  │   dplyr::left_join(horizon_data, sim_component_comp(component_data),
                  │                    by = "cokey")
                  ▼
    horizon-level data now carrying sim_comppct on every row
                  │
                  ▼

depth-simulation.R
======================
get_aws_data_by_mukey() ──► raw SSURGO horizon data (SDA_query)
                  │
                  ▼
   sim_component_comp() + dplyr::left_join(by = "cokey")   (performed internally
                  │                                         by simulate_profile_depths_by_mukey()
                  ▼                                         as of the sim_comppct integration fix)
   simulate_profile_depths_by_mukey()  ──┐
                  │                       │  (per-component profile loop)
                  ▼                       │
   simulate_profile_depths_by_collection()│
   simulate_profile_depths_by_collection_parallel()
                  │                       │
                  └───────────────────────┘
                  ▼
      simulate_and_perturb_soil_profiles()   *** still requires horizons$sim_comppct
                                                  to already be present - callers who use
                                                  this function directly (not via
                                                  simulate_profile_depths_by_mukey()) must
                                                  still perform the join themselves ***
                  │
                  ├──► simulate_soil_profile_thickness()
                  │        ├──► simulate_soil_profile_top_down()   [tri_dist()]
                  │        └──► simulate_soil_profile_bottom_up()  [tri_dist()]
                  │
                  ├──► query_osd_distinctness()
                  │        └──► infill_missing_distinctness()
                  │
                  ├──► infill_missing_depth_variability()   (via thickness step)
                  │
                  ├──► aqp::perturb() x2 (thickness, then boundary)
                  │
                  └──► adjust_out_of_range_profiles() [internal, clamps to l/h bounds]
                  ▼
      simulated SoilProfileCollection
                  │
                  ▼
      evaluate_simulated_depths()  ──► QC report of any still-out-of-range horizons
```

## Dependencies

- **`aqp`** - `SoilProfileCollection` construction/manipulation
  (`aqp::depths<-`, `aqp::horizons()`, `aqp::hzdesgnname<-`,
  `aqp::generalizeHz()`, `aqp::hzDistinctnessCodeToOffset()`,
  `aqp::perturb()`, `aqp::combine()`, `aqp::profile_id<-`) - used throughout
  `depth-simulation.R` only; `property-simulation.R` has no `aqp` dependency.
- **`soilDB`** - `soilDB::SDA_query()` (`get_aws_data_by_mukey()`) and
  `soilDB::fetchOSD()` (`query_osd_distinctness()`) for live SSURGO/OSD
  database access.
- **`future` / `future.apply`** - parallel dispatch backend for
  `simulate_profile_depths_by_collection_parallel()`
  (`future::plan(future::multisession)`, `future.apply::future_lapply()`).
- **`dplyr`** - data manipulation throughout both files (`group_by`,
  `mutate`, `left_join`, `select`, `summarise`, `rowwise`, etc.).
- **`R/distributions.R`** (soilSIM internal) - `tri_dist()` (used by both
  files: `sim_component_comp()`/`simulate_correlated_triangular()` in
  `property-simulation.R`; all four depth-simulation functions in
  `depth-simulation.R`) and `ilr_forward()`/`ilr_inverse()` (used only by
  `simulate_cokey_generalized()` for compositional texture handling).
- **`R/ssurgo-acquisition.R` / `R/ssurgo-processing.R`** (soilSIM internal) -
  upstream SSURGO data-source adapters; `get_aws_data_by_mukey()` in this
  group duplicates rather than reuses
  `execute_ssurgo_query_working()`'s query (see that function's docs for
  why), but both ultimately query the same underlying SSURGO tables.
- **Downstream consumers**:
  - **`R/ssurgo-simulation.R`'s `simulate_ssurgo_mapunit_draws()`** is a
    real, live in-package consumer of this group: it calls
    `remove_organic_layer()` and `sim_component_comp()` directly, performs
    the `sim_comppct` component-to-horizon join itself, then calls
    `simulate_cokey_generalized()` per cokey - i.e. it exercises exactly the
    `property-simulation.R` half of this group (component composition +
    per-cokey property simulation) as part of the multi-source raster
    fusion pipeline's SSURGO adapter (prior-side data for
    `R/raster-fusion.R`'s `fuse_property_adaptive()`/`fuse_texture_group()`).
    It does **not** call anything in `depth-simulation.R`.
  - **AWS / van Genuchten modeling** (`R/aws-simulation.R`) is documented in
    `R/soilSIM-package.R` as conceptually downstream in the pipeline (depth-
    sliced, depth-simulated profiles feeding available-water-storage
    estimation, per the legacy `code/soil_simulation_summary.md` workflow),
    but there is no direct function-level call between `depth-simulation.R`
    and `aws-simulation.R` in the current codebase - `aws-simulation.R`'s
    own header notes it is self-contained and does not depend on
    `sim_component_comp()` or any other helper from this group.

## Data Flow In/Out

**In:**
- SSURGO component data (`mukey`, `cokey`, `compname`, `comppct_l/r/h`) -
  input to `sim_component_comp()`.
- SSURGO horizon data (`hzname`, texture/bulk-density/water-retention/RFV/
  pH/CEC/OM `_l/_r/_h` triplets, `genhz`) - input to
  `simulate_cokey_generalized()`, `get_aws_data_by_mukey()`, and the depth-
  simulation functions.
- OSD (Official Series Description) horizon distinctness data, fetched
  live via `soilDB::fetchOSD()` inside `query_osd_distinctness()`.
- A `mukey` (string/numeric) - sole input to
  `simulate_profile_depths_by_mukey()`, which internally queries SSURGO.
- An existing `aqp::SoilProfileCollection` - input to
  `simulate_and_perturb_soil_profiles()` and the collection-level
  orchestration functions.

**Out:**
- Simulated component percentages: one row per component with a
  `sim_comppct` column (`sim_component_comp()`).
- Simulated per-cokey property realizations: one row per horizon x
  realization, with simulated bulk density/water retention/texture/RFV/pH/
  CEC/OM columns plus identifying columns
  (`simulate_cokey_generalized()`).
- Simulated horizon-thickness variability: one row per horizon with a
  `thickness_sd` estimate (`simulate_soil_profile_thickness()`).
- Simulated/perturbed `SoilProfileCollection`s: `n_simulations` realizations
  per input profile, with perturbed horizon top/bottom depths
  (`simulate_and_perturb_soil_profiles()` and its collection/mukey-level
  wrappers).
- QC data frames flagging out-of-range simulated depths
  (`evaluate_simulated_depths()`).

## Known Limitations

- **`simulate_and_perturb_soil_profiles()`** - requires `soil_profile`'s
  horizons to already carry a `sim_comppct` column (used to derive
  `n_simulations`); `sim_component_comp()` produces this column but at
  per-**component** grain, not per-**horizon** grain, so callers who use
  this function *directly* (not via `simulate_profile_depths_by_mukey()`)
  must still `dplyr::left_join()` the two by `cokey` before calling it, or
  it errors with a missing-column condition. This part of the contract is
  unchanged.
- ~~**`simulate_profile_depths_by_mukey()`** - has the same underlying
  `sim_comppct` requirement but does not itself call `sim_component_comp()`
  or perform the component-to-horizon join~~ **Fixed**: this function now
  calls `sim_component_comp(mu_data, n_simulations = n_simulations)` and
  left-joins the result onto every horizon row by `cokey` internally,
  mirroring the same pattern `R/ssurgo-simulation.R`'s
  `simulate_ssurgo_mapunit_draws()` already used for the property-simulation
  path. Its `n_simulations` parameter is no longer unused - it now flows
  directly into `sim_component_comp()`'s own `n_simulations` argument
  (number of triangular draws per component).

## Usage Example

```r
# simulate_profile_depths_by_mukey() now derives and joins sim_comppct
# internally (sim_component_comp() -> left_join by cokey), so a direct call
# is sufficient - no manual join step required:
simulated_profiles <- simulate_profile_depths_by_mukey(
  mukey = "123456", n_simulations = 1000, seed = 123
)

# Calling simulate_and_perturb_soil_profiles() directly (bypassing
# simulate_profile_depths_by_mukey()) still requires the manual join, since
# that lower-level function's own contract hasn't changed:
component_data <- sim_component_comp(ssurgo_component_data, n_simulations = 1000)
horizon_data <- dplyr::left_join(
  ssurgo_horizon_data,
  component_data[, c("cokey", "sim_comppct")],
  by = "cokey"
)
horizon_data$id <- horizon_data$compname
aqp::depths(horizon_data) <- id ~ hzdept_r + hzdepb_r
aqp::hzdesgnname(horizon_data) <- "hzname"

simulated_profiles <- simulate_and_perturb_soil_profiles(horizon_data[1, ])
```
