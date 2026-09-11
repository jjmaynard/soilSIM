# Simulate Soil Component Composition

Simulates component-composition percentages from low/rep/high
(`comppct_l/r/h`) values via
[`tri_dist()`](https://jjmaynard.github.io/soilSIM/reference/tri_dist.md),
one component (`cokey`) at a time, and derives a single `sim_comppct`
value per component.

## Usage

``` r
simulate_component_composition(data, n_simulations = 1000)
```

## Arguments

- data:

  Data frame with `mukey`, `cokey`, `compname`, `comppct_l`,
  `comppct_r`, `comppct_h` (one or more rows per component; deduplicated
  internally).

- n_simulations:

  Integer, number of triangular draws per component (default 1000).

## Value

A data frame, one row per component, with an added `sim_comppct` column.

## Grain mismatch with `R/core-simulation.R`

This returns one row per **component** (`cokey`), but
[`simulate_and_perturb_soil_profiles()`](https://jjmaynard.github.io/soilSIM/reference/simulate_and_perturb_soil_profiles.md)/
[`simulate_profile_depths()`](https://jjmaynard.github.io/soilSIM/reference/simulate_profile_depths.md)
need `sim_comppct` on every **horizon** row. Callers must
[`dplyr::left_join()`](https://dplyr.tidyverse.org/reference/mutate-joins.html)
this output onto horizon-level data by `cokey` before passing it to
those functions - this function alone does not satisfy their
requirement.

## `sim_comppct`'s derivation

`sim_comppct <- round(sum(<n_simulations> triangular draws of comppct) / 100)`,
which is `round(n_simulations * comppct_r / 100)` in expectation (e.g.
`comppct_r = 30`, `n_simulations = 1000` -\> `sim_comppct` ~= 300). It
is not an independently-meaningful simulation count.
