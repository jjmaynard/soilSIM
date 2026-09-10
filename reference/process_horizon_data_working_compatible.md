# Process Horizon Data (Working Compatible)

Process Horizon Data (Working Compatible)

## Usage

``` r
process_horizon_data_working_compatible(
  raw_data,
  detect_unsuitable = TRUE,
  advanced_cleaning = TRUE,
  standardize_names = TRUE,
  remove_invalid = TRUE,
  calculate_derived = TRUE,
  max_depth = DEFAULT_MAX_DEPTH_CM,
  verbose = FALSE
)
```

## Arguments

- raw_data:

  Raw SSURGO data containing horizon information

- detect_unsuitable:

  Logical; detect unsuitable horizons

- advanced_cleaning:

  Logical; apply advanced data cleaning functions

- standardize_names:

  Logical; standardize column names

- remove_invalid:

  Logical; remove invalid horizon records

- calculate_derived:

  Logical; calculate derived properties

- max_depth:

  Maximum depth for processing

- verbose:

  Logical; provide progress messages

## Value

List with processed horizon data and processing statistics
