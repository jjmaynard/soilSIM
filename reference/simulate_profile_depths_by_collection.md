# Run Soil Profile Depth Simulations for a SoilProfileCollection

Simulates soil profile depths for each profile in a
`SoilProfileCollection` by applying
[`simulate_and_perturb_soil_profiles()`](https://jjmaynard.github.io/soilSIM/reference/simulate_and_perturb_soil_profiles.md)
to each individual profile, then combines the results into a single
`SoilProfileCollection`.

## Usage

``` r
simulate_profile_depths_by_collection(soil_collection, seed = 123)
```

## Arguments

- soil_collection:

  A SoilProfileCollection object containing soil profile data.

- seed:

  An integer value for setting the random seed (default is 123) for
  reproducible simulations.

## Value

A SoilProfileCollection object containing the simulated soil profiles.

## Examples

``` r
if (FALSE) { # \dontrun{
  simulated_collection <- simulate_profile_depths_by_collection(soil_collection, seed = 123)
} # }
```
