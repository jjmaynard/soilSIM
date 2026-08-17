# Simulate Component Compositions (Enhanced)

Enhanced component composition simulation

## Usage

``` r
sim_component_compositions(
  component_data,
  n_realizations = 1000,
  config = NULL,
  verbose = getOption("ssurgo.verbose", FALSE)
)
```

## Arguments

- component_data:

  Data frame with component data

- n_realizations:

  Number of realizations

- config:

  Configuration settings

- verbose:

  Logical; if `TRUE`, temporarily raises the package's log level so
  `INFO`-level progress messages print for the duration of this call
  (default `FALSE` - quiet). See
  [`set_verbose_logging()`](https://jjmaynard.github.io/soilSIM/reference/set_verbose_logging.md).

## Value

List with simulated component compositions and quality metrics
