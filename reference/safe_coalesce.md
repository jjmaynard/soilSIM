# Safe Coalesce

Safely coalesces values from multiple columns, handling different data
types and missing values.

## Usage

``` r
safe_coalesce(data, column_names, target_type = "numeric", default_value = NA)
```

## Arguments

- data:

  Input data

- column_names:

  Vector of column names in priority order

- target_type:

  Target data type ("numeric", "character", "logical")

- default_value:

  Default value if all columns are NA

## Value

Coalesced values
