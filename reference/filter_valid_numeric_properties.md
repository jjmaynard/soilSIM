# Filter Valid Numeric Properties

FIXED: Remove empty or all-NA columns to prevent min/max warnings

## Usage

``` r
filter_valid_numeric_properties(data, properties, verbose = FALSE)
```

## Arguments

- data:

  Input data

- properties:

  Properties to check

- verbose:

  Logical; provide progress messages

## Value

Filtered properties list
