# Bottom-up Simulation of Soil Profile

Simulates a soil profile from the bottom upward by calculating the top
and bottom depths for each horizon. The simulation starts from the
bottom horizon and proceeds upward, setting the bottom depth of each
horizon to the top depth of the horizon immediately below. For each
horizon, the bottom and top depths are simulated using a triangular
distribution via
[`tri_dist()`](https://jjmaynard.github.io/soilSIM/reference/tri_dist.md).
The function ensures that there are no gaps between horizons and that
the top of the uppermost horizon is set to 0.

## Usage

``` r
simulate_soil_profile_bottom_up(horizon_data)
```

## Arguments

- horizon_data:

  A data frame containing horizon data. The data frame must include the
  following columns:

  - `hzname`: The name of the horizon.

  - `hzdept_l`: The lower bound for the top depth of the horizon.

  - `hzdept_r`: The representative top depth of the horizon.

  - `hzdept_h`: The upper bound for the top depth of the horizon.

  - `hzdepb_l`: The lower bound for the bottom depth of the horizon.

  - `hzdepb_r`: The representative bottom depth of the horizon.

  - `hzdepb_h`: The upper bound for the bottom depth of the horizon.

## Value

A data frame with the following columns:

- `hzname`: The name of the horizon.

- `top`: The simulated top depth of the horizon.

- `bottom`: The simulated bottom depth of the horizon.

## Examples

``` r
horizon_data <- data.frame(
  hzname = c("A", "B", "C"),
  hzdept_l = c(0, 20, 40),
  hzdept_r = c(0, 25, 45),
  hzdept_h = c(0, 30, 50),
  hzdepb_l = c(20, 40, 60),
  hzdepb_r = c(25, 45, 65),
  hzdepb_h = c(30, 50, 70)
)
profile <- simulate_soil_profile_bottom_up(horizon_data)
print(profile)
#>   hzname top bottom
#> 1      A   0     23
#> 2      B  23     47
#> 3      C  47     62
```
