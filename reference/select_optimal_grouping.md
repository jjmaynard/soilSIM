# Select Optimal Grouping Strategy

Enhanced version that uses Module 0 validation utilities

## Usage

``` r
select_optimal_grouping(
  data,
  min_profiles,
  min_obs,
  target_groups,
  verbose = getOption("ssurgo.verbose", FALSE)
)
```

## Arguments

- data:

  Input NRCS data

- min_profiles:

  Minimum profiles per group

- min_obs:

  Minimum observations per group

- target_groups:

  Minimum number of adequate groups desired

- verbose:

  Logical; if `TRUE`, temporarily raises the package's log level so
  `INFO`-level progress messages print for the duration of this call
  (default `FALSE` - quiet). See
  [`set_verbose_logging()`](https://jjmaynard.github.io/soilSIM/reference/set_verbose_logging.md).

## Value

Selected grouping strategy name
