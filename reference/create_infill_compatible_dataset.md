# Create Infill-Compatible Dataset

Creates main processed dataset that's fully compatible with
infill_soil_property() Uses Module 8 utilities for data manipulation

## Usage

``` r
create_infill_compatible_dataset(
  raw_data,
  horizon_processing,
  component_processing,
  options,
  verbose = FALSE
)
```

## Arguments

- raw_data:

  Original raw data

- horizon_processing:

  Horizon processing results

- component_processing:

  Component processing results

- options:

  Processing options

- verbose:

  Logical; provide progress messages

## Value

Data frame compatible with existing infill workflow
