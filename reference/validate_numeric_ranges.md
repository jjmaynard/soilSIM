# Validate Numeric Ranges

Validates that numeric properties fall within realistic ranges.

## Usage

``` r
validate_numeric_ranges(data, property_ranges, action = "warn")
```

## Arguments

- data:

  Input data

- property_ranges:

  Named list of property ranges

- action:

  Action for out-of-range values: "warn", "error", or "clip"

## Value

Range validation results
