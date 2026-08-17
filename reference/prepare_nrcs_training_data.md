# Prepare NRCS Training Data for GP Model Building

Prepare NRCS Training Data for GP Model Building

## Usage

``` r
prepare_nrcs_training_data(
  nrcs_combined_data,
  grouping_strategy = "auto",
  min_profiles_per_group = 3,
  min_observations_per_group = 15,
  target_min_groups = 3,
  max_depth = 250,
  validation_config = NULL,
  verbose = getOption("ssurgo.verbose", FALSE)
)
```

## Arguments

- nrcs_combined_data:

  Combined NRCS data with horizon and component information

- grouping_strategy:

  How to group soils ("auto" for intelligent selection)

- min_profiles_per_group:

  Minimum profiles required per group (default = 3)

- min_observations_per_group:

  Minimum total observations per group (default = 15)

- target_min_groups:

  Minimum number of adequate groups desired (default = 3)

- max_depth:

  Maximum depth to include in training (default = 250 cm)

- validation_config:

  Validation configuration from Module 0

- verbose:

  Logical; if `TRUE`, temporarily raises the package's log level so
  `INFO`-level progress messages print for the duration of this call
  (default `FALSE` - quiet). See
  [`set_verbose_logging()`](https://jjmaynard.github.io/soilSIM/reference/set_verbose_logging.md).

## Value

Processed data frame ready for GP model building
