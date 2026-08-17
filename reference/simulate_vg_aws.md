# Simulate Available Water Holding Capacity from van Genuchten Parameters

Monte Carlo-simulates field capacity (FC) and permanent wilting point
(PWP) water content - and their difference, available water holding
capacity (AWHC) - for each row (soil layer/component) of van Genuchten
parameters, by drawing `n_simulations` random samples per parameter from
`stats::rnorm(mean, sd)` and evaluating
[`van_genuchten()`](https://jjmaynard.github.io/soilSIM/reference/van_genuchten.md)
at the FC/PWP matric potentials.

## Usage

``` r
simulate_vg_aws(data, n_simulations = 100)
```

## Arguments

- data:

  A data frame with one row per soil layer/component, and columns
  `alpha`, `sd_alpha`, `npar`, `sd_npar`, `theta_r`, `sd_theta_r`,
  `theta_s`, `sd_theta_s`, `layerID` - exactly the shape produced by
  `soilDB::ROSETTA(..., include.sd = TRUE)` plus a caller-added
  `layerID` (see
  [`calculate_aws_df()`](https://jjmaynard.github.io/soilSIM/reference/calculate_aws_df.md)).

- n_simulations:

  Number of Monte Carlo draws per row (default 100).

## Value

A named list (one element per row, keyed
`paste(layerID, i, sep = "_")`), each a data frame of `n_simulations`
draws with columns `alpha`, `n`, `theta_r`, `theta_s`, `sim_num`,
`theta_fc`, `theta_pwp`, `AWHC`. Rows with any missing van Genuchten
parameter are skipped (absent from the result).

## Ported-as-is behavioral quirks

Two aspects of this function are preserved exactly as in the legacy
source rather than "fixed," since both are plausibly intentional design
choices (not unconditional-crash bugs) and changing either would alter
the actual simulated values every caller receives:

- `set.seed(123)` is called **inside** the per-row loop, so every row
  re-seeds and draws from an identically-seeded random stream rather
  than an evolving one across rows.

- Sampled `alpha`/`n` values are back-transformed via `10^(...)`, i.e.
  the input `data$alpha`/`data$npar` (and their SDs) are treated as
  already being in log10 space. This is a standard technique for keeping
  van Genuchten shape parameters positive after adding Gaussian noise,
  but whether
  [`soilDB::ROSETTA()`](http://ncss-tech.github.io/soilDB/reference/ROSETTA.md)'s
  actual reported values are meant to be interpreted this way is not
  verified by this port.
