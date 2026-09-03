# Check Required Columns

Validates that required columns exist and have appropriate data types.

## Usage

``` r
check_required_columns(data, column_specifications, strict_mode = TRUE)
```

## Arguments

- data:

  Input data

- column_specifications:

  List of column specifications

- strict_mode:

  Whether to enforce strict type checking

## Value

Column validation results

## Usage note

Fully implemented and exported, but not currently called from anywhere
else in the package (confirmed by source grep) - standalone public API
for callers who need it, not dead/broken code.
