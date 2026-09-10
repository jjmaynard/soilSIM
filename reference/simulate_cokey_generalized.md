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
  txt_correlation_matrices = NULL,
  requested_properties = NULL
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
  `c("db", "wr_3b", "wr_15b", "ilr1", "ilr2", "rfv", "ph", "cec", "soc", "caco3", "ec", "ecec", "gypsum", "sar")`.
  The 5 chemistry properties `caco3`/`ec`/`ecec`/`gypsum`/`sar` have no
  real KSSL-fit correlation data (see
  [`build_kssl_fallback_matrix()`](https://jjmaynard.github.io/soilSIM/reference/build_kssl_fallback_matrix.md)) -
  [`simulate_ssurgo_mapunit_draws()`](https://jjmaynard.github.io/soilSIM/reference/simulate_ssurgo_mapunit_draws.md)
  builds its `correlation_matrices` through that function so they're
  present (as identity/uncorrelated) rather than missing; a direct
  caller supplying its own matrix without them will simply never
  populate `param_list[["caco3"]]` etc. (see `get_param_set()`'s per-row
  gating below), not error.

- txt_correlation_matrices:

  A list of texture correlation matrices keyed by `genhz`
  (sand/silt/clay order, matching the
  [`simulate_correlated_triangular()`](https://jjmaynard.github.io/soilSIM/reference/simulate_correlated_triangular.md)
  call for texture).

- requested_properties:

  `NULL` (default) simulates every property with a complete `_l/_r/_h`
  triplet, unchanged. Otherwise a character vector of property
  identifiers - in any vocabulary
  [`normalize_requested_properties()`](https://jjmaynard.github.io/soilSIM/reference/normalize_requested_properties.md)
  accepts (caller ids, SSURGO stems, output column names, or the
  internal `param_order` names) - restricting the simulation to just
  those (plus the texture coupling: any texture member pulls in the full
  sand/silt/clay draw and both ILR axes). The Gaussian-copula marginals
  of a restricted simulation match the corresponding columns of a full
  one under the default `joint_copula` vertical-correlation method
  (modulo RNG-stream Monte Carlo noise); see
  [`simulate_ssurgo_mapunit_draws()`](https://jjmaynard.github.io/soilSIM/reference/simulate_ssurgo_mapunit_draws.md)
  for the `gp_quantile_retrofit` caveat.

## Value

A data frame of simulated property values across all rows/realizations,
with `compname`, `mukey`, `cokey`, `hzdept_r`, `hzdepb_r`,
`simulation_number`, `unique_id`, and (when `sim_cokey` itself has a
`bound_sd` column - see
[`attach_osd_boundary_distinctness()`](https://jjmaynard.github.io/soilSIM/reference/attach_osd_boundary_distinctness.md))
`bound_sd`. The `soc` column is an SSURGO-organic-matter-**derived SOC
estimate** (`om * OM_TO_SOC_FACTOR`, the inverse Van Bemmelen factor) -
not a lab-measured soil organic carbon value. See `OM_TO_SOC_FACTOR`'s
own docs.
