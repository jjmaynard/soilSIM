# Export Workflow Metadata

Exports comprehensive metadata about the workflow execution.

## Usage

``` r
export_workflow_metadata(
  workflow_results,
  output_path,
  include_data_summary = TRUE,
  verbose = getOption("ssurgo.verbose", FALSE)
)
```

## Arguments

- workflow_results:

  Workflow results

- output_path:

  Output file path

- include_data_summary:

  Whether to include data summaries

- verbose:

  Logical; if `TRUE`, temporarily raises the package's log level so
  `INFO`-level progress messages print for the duration of this call
  (default `FALSE` - quiet). See
  [`set_verbose_logging()`](https://jjmaynard.github.io/soilSIM/reference/set_verbose_logging.md).

## Value

Success status
