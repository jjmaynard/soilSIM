# Convert to Property Matrices

Enhanced version with Module 8 safe operations and validation.

## Usage

``` r
convert_to_property_matrices(
  simulation_data,
  properties,
  unique_depths,
  sim_numbers,
  verbose = getOption("ssurgo.verbose", FALSE)
)
```

## Arguments

- simulation_data:

  Simulation data in long format

- properties:

  Properties to convert

- unique_depths:

  Depth vector

- sim_numbers:

  Simulation numbers

- verbose:

  Logical; if `TRUE`, temporarily raises the package's log level so
  `INFO`-level progress messages print for the duration of this call
  (default `FALSE` - quiet). See
  [`set_verbose_logging()`](https://jjmaynard.github.io/soilSIM/reference/set_verbose_logging.md).

## Value

Named list of property matrices
