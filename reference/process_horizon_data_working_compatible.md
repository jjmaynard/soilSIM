# Process Horizon Data (Working Compatible)

Enhanced horizon processing that uses Module 8 utilities and proven
working functions

## Usage

``` r
process_horizon_data_working_compatible(
  raw_data,
  detect_unsuitable = TRUE,
  advanced_cleaning = TRUE,
  standardize_names = TRUE,
  remove_invalid = TRUE,
  calculate_derived = TRUE,
  max_depth = 250,
  verbose = FALSE
)
```

## Arguments

- raw_data:

  Raw SSURGO data containing horizon information

- detect_unsuitable:

  Logical; detect unsuitable horizons using Module 8 function

- advanced_cleaning:

  Logical; apply advanced data cleaning functions

- standardize_names:

  Logical; standardize column names using Module 8

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
