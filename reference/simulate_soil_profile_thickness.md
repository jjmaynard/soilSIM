# Simulate Soil Profile Thickness

Simulates the thickness of soil horizons for a given soil profile using
both top-down and bottom-up simulation methods. It performs multiple
simulations and calculates the standard deviation of horizon thickness
(i.e., the difference between representative bottom and top depths)
across the simulations.

## Usage

``` r
simulate_soil_profile_thickness(horizon_data, n_simulations = 500)
```

## Arguments

- horizon_data:

  A dataframe containing horizon data. It must include the following
  columns: - `hzname`: The name of each horizon (e.g., A, B, C). -
  `hzdept_r`: The representative top depth of each horizon. -
  `hzdepb_r`: The representative bottom depth of each horizon. -
  `hzdept_l`: The low estimate for the top depth of each horizon. -
  `hzdept_h`: The high estimate for the top depth of each horizon. -
  `hzdepb_l`: The low estimate for the bottom depth of each horizon. -
  `hzdepb_h`: The high estimate for the bottom depth of each horizon.

- n_simulations:

  Integer, number of simulations to perform (default is 500).

## Value

A dataframe containing the following columns: - `hzname`: The horizon
name. - `top`: The representative top depth of the horizon (from
`hzdept_r`). - `bottom`: The representative bottom depth of the horizon
(from `hzdepb_r`). - `thickness_sd`: The standard deviation of the
horizon thickness (bottom - top) across the simulations.

## Examples

``` r
horizon_data <- data.frame(
  hzname = c("A", "B", "C"),
  hzdept_r = c(0, 20, 35),
  hzdepb_r = c(20, 35, 50),
  hzdept_l = c(0, 15, 30),
  hzdept_h = c(0, 25, 40),
  hzdepb_l = c(15, 30, 45),
  hzdepb_h = c(25, 40, 60)
)
set.seed(123)
summarized_results <- simulate_soil_profile_thickness(horizon_data, n_simulations = 500)
print(summarized_results)
#>   hzname top bottom thickness_sd
#> 1      A   0     20         2.00
#> 2      B  20     35         2.91
#> 3      C  35     50         3.69
```
