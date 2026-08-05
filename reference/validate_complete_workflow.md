# Validate Complete Workflow

Master validation function enhanced for comprehensive assessment of the
entire soil simulation workflow.

## Usage

``` r
validate_complete_workflow(
  workflow_results,
  original_data = NULL,
  validation_config = NULL,
  generate_plots = TRUE,
  output_dir = NULL,
  verbose = getOption("ssurgo.verbose", FALSE)
)
```

## Arguments

- workflow_results:

  Complete workflow results including all intermediate steps

- original_data:

  Original SSURGO/NRCS data for comparison

- validation_config:

  Configuration for validation parameters (uses Module 0 defaults)

- generate_plots:

  Whether to generate diagnostic plots (default = TRUE)

- output_dir:

  Directory for saving validation outputs

- verbose:

  Logical; if `TRUE`, temporarily raises the package's log level so
  `INFO`-level progress messages print for the duration of this call
  (default `FALSE` - quiet). See
  [`set_verbose_logging()`](https://jjmaynard.github.io/soilSIM/reference/set_verbose_logging.md).

## Value

Comprehensive validation results
