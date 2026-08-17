# Correct Distribution Shapes

Enhanced version with Module 8 property validation and constraints.

## Usage

``` r
correct_distribution_shapes(
  adjusted_data,
  original_data,
  properties,
  config = NULL,
  verbose = getOption("ssurgo.verbose", FALSE)
)
```

## Arguments

- adjusted_data:

  Adjusted simulation data

- original_data:

  Original simulation data

- properties:

  Properties to validate

- config:

  Configuration settings

- verbose:

  Logical; if `TRUE`, temporarily raises the package's log level so
  `INFO`-level progress messages print for the duration of this call
  (default `FALSE` - quiet). See
  [`set_verbose_logging()`](https://jjmaynard.github.io/soilSIM/reference/set_verbose_logging.md).

## Value

Corrected simulation data
