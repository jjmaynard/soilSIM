# Validate Correlation Preservation

Helper function from GP-depth-adjust.R to validate correlation
preservation

## Usage

``` r
validate_correlation_preservation(
  original_list,
  adjusted_list,
  depths,
  verbose = getOption("ssurgo.verbose", FALSE)
)
```

## Arguments

- original_list:

  Original simulated property list

- adjusted_list:

  Adjusted property list

- depths:

  Depth vector

- verbose:

  Logical; if `TRUE`, temporarily raises the package's log level so
  `INFO`-level progress messages print for the duration of this call
  (default `FALSE` - quiet). See
  [`set_verbose_logging()`](https://jjmaynard.github.io/soilSIM/reference/set_verbose_logging.md).
