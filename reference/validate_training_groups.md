# Validate Training Groups

Enhanced validation using Module 0 utilities

## Usage

``` r
validate_training_groups(
  processed_data,
  properties = c("clay_pct", "sand_pct", "pH", "organic_matter"),
  min_profiles = 3,
  min_observations = 15,
  verbose = getOption("ssurgo.verbose", FALSE)
)
```

## Arguments

- processed_data:

  Processed NRCS data

- properties:

  Properties to validate

- min_profiles:

  Minimum profiles per group

- min_observations:

  Minimum observations per group

- verbose:

  Logical; if `TRUE`, temporarily raises the package's log level so
  `INFO`-level progress messages print for the duration of this call
  (default `FALSE` - quiet). See
  [`set_verbose_logging()`](https://jjmaynard.github.io/soilSIM/reference/set_verbose_logging.md).

## Value

Validation results list
