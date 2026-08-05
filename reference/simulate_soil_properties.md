# Enhanced Soil Property Simulation with NRCS GP Models + Cholesky Correlations

Complete implementation from GP-depth-adjust.R that combines triangular
distribution

- Cholesky approach with NRCS GP models for realistic depth trends

## Usage

``` r
simulate_soil_properties(
  target_cokey,
  nrcs_gp_models = NULL,
  cokey_mapping = NULL,
  sim_data,
  correlation_matrices,
  txt_correlation_matrices,
  n_simulations = 100,
  use_nrcs_gp = !is.null(nrcs_gp_models),
  preserve_correlations = TRUE,
  verbose = getOption("ssurgo.verbose", FALSE)
)
```

## Arguments

- target_cokey:

  Character string of the cokey to simulate

- nrcs_gp_models:

  Your fitted NRCS GP models (optional)

- cokey_mapping:

  Mapping from match_simulated_soils_to_gp_models (optional)

- sim_data:

  Your processed simulation data for this cokey

- correlation_matrices:

  Your existing correlation matrices by genetic horizon

- txt_correlation_matrices:

  Your existing texture correlation matrices

- n_simulations:

  Number of Monte Carlo realizations

- use_nrcs_gp:

  Whether to use NRCS GP models (TRUE) or local GP models (FALSE)

- preserve_correlations:

  Whether to preserve within-depth correlations

- verbose:

  Logical; if `TRUE`, temporarily raises the package's log level so
  `INFO`-level progress messages print for the duration of this call
  (default `FALSE` - quiet). See
  [`set_verbose_logging()`](https://jjmaynard.github.io/soilSIM/reference/set_verbose_logging.md).

## Value

Enhanced simulation results with realistic depth trends and preserved
correlations
