# Top-down Simulation of Soil Profile

Simulates a soil profile by calculating the top and bottom depths for
each horizon in a top-down manner. The first horizon always starts at 0
cm, and the top of each subsequent horizon is set equal to the bottom
depth of the previous horizon. The bottom depth of each horizon is
simulated using a triangular distribution via
[`tri_dist()`](https://jjmaynard.github.io/soilSIM/reference/tri_dist.md),
defined by a lower bound (`hzdepb_l`), a representative value
(`hzdepb_r`), and an upper bound (`hzdepb_h`).

## Usage

``` r
simulate_soil_profile_top_down(horizon_data)
```

## Arguments

- horizon_data:

  A data frame containing horizon data. It should include the following
  columns:

  - `hzname`: Name of the horizon.

  - `hzdepb_l`: Lower bound for the bottom depth of the horizon.

  - `hzdepb_r`: Representative bottom depth of the horizon.

  - `hzdepb_h`: Upper bound for the bottom depth of the horizon.

## Value

A data frame with the following columns:

- `hzname`: The name of the horizon.

- `top`: The simulated top depth of the horizon.

- `bottom`: The simulated bottom depth of the horizon.
