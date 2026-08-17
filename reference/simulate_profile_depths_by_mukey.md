# Simulate Profile Depths by Mukey

Queries soil data for a given map unit key (mukey) using the SSURGO
database, sets up the corresponding soil profile data, and then
simulates and perturbs the soil profiles for that mukey via
[`simulate_and_perturb_soil_profiles()`](https://jjmaynard.github.io/soilSIM/reference/simulate_and_perturb_soil_profiles.md).

## Usage

``` r
simulate_profile_depths_by_mukey(mukey, n_simulations = 100, seed = 123)
```

## Arguments

- mukey:

  A string or numeric value representing the map unit key to query.

- n_simulations:

  Integer, the number of triangular draws per component used to derive
  `sim_comppct` via
  [`sim_component_comp()`](https://jjmaynard.github.io/soilSIM/reference/sim_component_comp.md)
  (default is 100).

- seed:

  An integer to set the random seed for reproducibility (default is
  123).

## Value

A SoilProfileCollection object containing the simulated and perturbed
soil profiles.

## Examples

``` r
if (FALSE) { # \dontrun{
  simulated_profiles <- simulate_profile_depths_by_mukey("123456", n_simulations = 100, seed = 123)
  print(simulated_profiles)
} # }
```
