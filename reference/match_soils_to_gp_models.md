# Match Soils to GP Models

Enhanced version with better error handling and fallback strategies

## Usage

``` r
match_soils_to_gp_models(
  simulated_cokeys,
  nrcs_combined_data,
  gp_models,
  property = "clay_pct",
  matching_strategy = "exact_cokey",
  verbose = getOption("ssurgo.verbose", FALSE)
)
```

## Arguments

- simulated_cokeys:

  Vector of cokeys from simulation data

- nrcs_combined_data:

  Original NRCS data used for GP training

- gp_models:

  Fitted GP models from build_stratified_gp_models()

- property:

  Property name for matching (default = "clay_pct")

- matching_strategy:

  Matching approach (default = "exact_cokey")

- verbose:

  Logical; if `TRUE`, temporarily raises the package's log level so
  `INFO`-level progress messages print for the duration of this call
  (default `FALSE` - quiet). See
  [`set_verbose_logging()`](https://jjmaynard.github.io/soilSIM/reference/set_verbose_logging.md).

## Value

Data frame mapping cokeys to GP model groups
