# Simulate Soil Properties for a Specific Cokey (Generalized)

Simulates a flexible set of soil properties (whichever of bulk density,
water retention, texture, rock fragment volume, pH, CEC, organic carbon
have complete `_l/_r/_h` triplets present) for each row of a single
cokey's horizon data, via
[`simulate_correlated_triangular()`](https://jjmaynard.github.io/soilSIM/reference/simulate_correlated_triangular.md),
using a genhz-keyed local correlation matrix. If texture columns are
present, sand/silt/clay are simulated jointly and converted to ILR
coordinates (`ilr1`/`ilr2`) before being included in the main correlated
simulation, then converted back after.

## Usage

``` r
simulate_cokey_generalized(
  sim_cokey,
  correlation_matrices,
  txt_correlation_matrices = NULL
)
```

## Arguments

- sim_cokey:

  A data frame of horizon rows for one cokey, with a `genhz` column and
  a `sim_comppct` column (number of realizations to simulate for that
  row - see
  [`sim_component_comp()`](https://jjmaynard.github.io/soilSIM/reference/sim_component_comp.md)).

- correlation_matrices:

  A list of correlation matrices keyed by `genhz`, with row/column names
  matching (a subset of)
  `c("db", "wr_3b", "wr_15b", "ilr1", "ilr2", "rfv", "ph", "cec", "soc")`.

- txt_correlation_matrices:

  A list of texture correlation matrices keyed by `genhz`
  (sand/silt/clay order, matching the
  [`simulate_correlated_triangular()`](https://jjmaynard.github.io/soilSIM/reference/simulate_correlated_triangular.md)
  call for texture).

## Value

A data frame of simulated property values across all rows/realizations,
with `compname`, `mukey`, `cokey`, `hzdept_r`, `hzdepb_r`,
`simulation_number`, `unique_id`.
