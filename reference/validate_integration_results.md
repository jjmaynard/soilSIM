# Validate Integration Results

Enhanced validation using Module 8 validation framework.

## Usage

``` r
validate_integration_results(
  original_data,
  integrated_data,
  properties,
  preserve_correlations = TRUE,
  verbose = getOption("ssurgo.verbose", FALSE)
)
```

## Arguments

- original_data:

  Original simulation data

- integrated_data:

  Integrated simulation data

- properties:

  Properties that were processed

- preserve_correlations:

  Whether correlations should be preserved

- verbose:

  Logical; if `TRUE`, temporarily raises the package's log level so
  `INFO`-level progress messages print for the duration of this call
  (default `FALSE` - quiet). See
  [`set_verbose_logging()`](https://jjmaynard.github.io/soilSIM/reference/set_verbose_logging.md).

## Value

Validation results list
