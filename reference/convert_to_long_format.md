# Convert to Long Format

Reshapes a named list of adjusted depth-by-simulation property matrices
back into long format, re-attaching component metadata.

## Usage

``` r
convert_to_long_format(
  adjusted_matrices,
  unique_depths,
  sim_numbers,
  original_data,
  properties,
  verbose = getOption("ssurgo.verbose", FALSE)
)
```

## Arguments

- adjusted_matrices:

  List of adjusted property matrices

- unique_depths:

  Depth vector

- sim_numbers:

  Simulation numbers

- original_data:

  Original simulation data for metadata

- properties:

  Properties that were adjusted

- verbose:

  Logical; if `TRUE`, temporarily raises the package's log level so
  `INFO`-level progress messages print for the duration of this call
  (default `FALSE` - quiet). See
  [`set_verbose_logging()`](https://jjmaynard.github.io/soilSIM/reference/set_verbose_logging.md).

## Value

Data frame in long format
