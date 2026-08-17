# Integrate Monte Carlo Simulations with GP Models

Master function that integrates Monte Carlo simulation results with GP
models to apply realistic depth trends while preserving within-depth
correlations. Enhanced with Module 8 utilities for robust processing and
validation.

## Usage

``` r
integrate_monte_carlo_with_gp(
  simulation_results,
  gp_models = NULL,
  cokey_mapping = NULL,
  integration_method = "hybrid",
  preserve_correlations = TRUE,
  properties = NULL,
  parallel = FALSE,
  n_cores = NULL,
  config = NULL,
  verbose = getOption("ssurgo.verbose", FALSE)
)
```

## Arguments

- simulation_results:

  Results from monte_carlo::generate_monte_carlo_realizations()

- gp_models:

  Optional NRCS GP models from gp_modeling module

- cokey_mapping:

  Optional mapping from simulations to NRCS GP groups

- integration_method:

  Method: "nrcs_gp", "local_gp", or "hybrid" (default = "hybrid")

- preserve_correlations:

  Whether to preserve within-depth correlations (default = TRUE)

- properties:

  Properties to adjust (NULL = auto-detect)

- parallel:

  Whether to use parallel processing (default = FALSE)

- n_cores:

  Number of cores for parallel processing

- config:

  Integration configuration (uses Module 8 defaults if NULL)

- verbose:

  Logical; if `TRUE`, temporarily raises the package's log level so
  `INFO`-level progress messages print for the duration of this call
  (default `FALSE` - quiet). See
  [`set_verbose_logging()`](https://jjmaynard.github.io/soilSIM/reference/set_verbose_logging.md).

## Value

Integrated simulation results with realistic depth trends
