# Generate Validation Report

Enhanced report generation I/O utilities and error handling.

## Usage

``` r
generate_validation_report(
  validation_results,
  output_format = "html",
  output_file = NULL,
  include_plots = TRUE,
  verbose = getOption("ssurgo.verbose", FALSE)
)
```

## Arguments

- validation_results:

  Results from validate_complete_workflow()

- output_format:

  Format for report: "html", "pdf", or "markdown"

- output_file:

  Output file path

- include_plots:

  Whether to include diagnostic plots

- verbose:

  Logical; if `TRUE`, temporarily raises the package's log level so
  `INFO`-level progress messages print for the duration of this call
  (default `FALSE` - quiet). See
  [`set_verbose_logging()`](https://jjmaynard.github.io/soilSIM/reference/set_verbose_logging.md).

## Value

Success status
