# Simulate and Perturb Soil Profiles

Takes a `SoilProfileCollection` object and perturbs its horizons to
simulate variability in soil profile thickness and boundary depths. The
simulation is performed in two steps: first, horizon thickness is
perturbed based on simulated thickness variability (using a top-down or
bottom-up simulation method), and then the boundary depths are perturbed
using distinctness-derived offsets. If the profile contains only one
horizon, the function simply replicates the profile.

## Usage

``` r
simulate_and_perturb_soil_profiles(soil_profile)
```

## Arguments

- soil_profile:

  A SoilProfileCollection object containing soil profile data.

## Value

A SoilProfileCollection object with perturbed horizon depths.

## Known limitation

This function requires `soil_profile`'s horizons to already carry a
`sim_comppct` column (it derives the number of simulations to run from
`unique(horizons(soil_profile)$sim_comppct)`).
[`sim_component_comp()`](https://jjmaynard.github.io/soilSIM/reference/sim_component_comp.md)
(`R/property-simulation.R`) produces this column, but at component
(`cokey`) grain, not horizon grain - callers must
`dplyr::left_join(horizon_data, sim_component_comp(component_data), by = "cokey")`
before calling this function; it will still error with a missing-column
condition if that join hasn't been done first.
