# Process SSURGO Data (Main Entry Point)

Enhanced version that maintains compatibility with working
infill_soil_property() while adding comprehensive processing,
validation, and reporting capabilities. Uses Module 8 utilities for
general processing tasks.

## Usage

``` r
process_ssurgo_data(
  raw_data,
  processing_options = list(),
  validate_results = TRUE,
  max_depth = 250,
  verbose = getOption("ssurgo.verbose", FALSE)
)
```

## Arguments

- raw_data:

  Raw SSURGO data from download functions

- processing_options:

  List of processing options and parameters

- validate_results:

  Logical; validate processing results (default: TRUE)

- max_depth:

  Maximum depth (cm) for processing (default: 250)

- verbose:

  Logical; provide detailed progress messages

## Value

List containing:

- processed_data: Combined processed data (compatible with infill
  functions)

- horizon_data: Processed horizon data

- component_data: Processed component data

- processing_metadata: Processing information and statistics

- validation_results: Data validation results

- quality_report: Comprehensive quality assessment
