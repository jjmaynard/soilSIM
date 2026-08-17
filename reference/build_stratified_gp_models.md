# Build Stratified GP Models

Enhanced version with proper Module 0 error handling and validation

## Usage

``` r
build_stratified_gp_models(
  processed_nrcs_data,
  properties = c("clay_pct", "sand_pct", "pH", "organic_matter"),
  min_profiles_per_group = 3,
  min_observations_per_group = 15,
  optimize_hyperparameters = TRUE,
  verbose = getOption("ssurgo.verbose", FALSE)
)
```

## Arguments

- processed_nrcs_data:

  Processed NRCS data from prepare_nrcs_training_data()

- properties:

  Character vector of properties to model

- min_profiles_per_group:

  Minimum profiles per group (default = 3)

- min_observations_per_group:

  Minimum observations per group (default = 15)

- optimize_hyperparameters:

  Whether to optimize GP hyperparameters (default = TRUE)

- verbose:

  Logical; if `TRUE`, temporarily raises the package's log level so
  `INFO`-level progress messages print for the duration of this call
  (default `FALSE` - quiet). See
  [`set_verbose_logging()`](https://jjmaynard.github.io/soilSIM/reference/set_verbose_logging.md).

## Value

List of GP models organized by property and group
