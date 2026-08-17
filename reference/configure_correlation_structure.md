# Configure Correlation Structure (Enhanced)

Enhanced correlation configuration

## Usage

``` r
configure_correlation_structure(
  simulation_params,
  properties,
  correlation_matrix = NULL,
  config,
  simulation_data = NULL,
  verbose = getOption("ssurgo.verbose", FALSE)
)
```

## Arguments

- simulation_params:

  Simulation parameters

- properties:

  Property names

- correlation_matrix:

  Optional correlation matrix

- config:

  Configuration settings

- simulation_data:

  Optional raw per-horizon data frame (passed through to
  [`estimate_property_correlations()`](https://jjmaynard.github.io/soilSIM/reference/estimate_property_correlations.md)
  so genhz-stratified estimation can be used when a `genhz` column is
  present).

- verbose:

  Logical; if `TRUE`, temporarily raises the package's log level so
  `INFO`-level progress messages print for the duration of this call
  (default `FALSE` - quiet). See
  [`set_verbose_logging()`](https://jjmaynard.github.io/soilSIM/reference/set_verbose_logging.md).

## Value

Enhanced correlation configuration
