# Calculate Workflow Performance

Derives real performance metrics from `components` (fraction of expected
components that were successfully extracted, i.e. non-`NULL`) and
`validation_results` (fraction of `*_validation` steps that neither were
skipped nor failed).

## Usage

``` r
calculate_workflow_performance(components, validation_results)
```

## Arguments

- components:

  Extracted workflow components (named list, possibly with `NULL`
  entries for components that could not be extracted).

- validation_results:

  Accumulated validation results, with one `<step>_validation` entry per
  pipeline step.

## Value

List with `processing_efficiency`, `validation_coverage`,
`overall_performance`.
