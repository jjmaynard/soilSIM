# Setup Distributions (Enhanced)

Enhanced distribution setup

## Usage

``` r
setup_distributions(
  simulation_params,
  properties,
  config,
  verbose = getOption("ssurgo.verbose", FALSE)
)
```

## Arguments

- simulation_params:

  List of simulation parameters by horizon

- properties:

  Character vector of property names

- config:

  Configuration settings

- verbose:

  Logical; if `TRUE`, temporarily raises the package's log level so
  `INFO`-level progress messages print for the duration of this call
  (default `FALSE` - quiet). See
  [`set_verbose_logging()`](https://jjmaynard.github.io/soilSIM/reference/set_verbose_logging.md).

## Value

Enhanced distribution configuration
