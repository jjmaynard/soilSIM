# Apply Validation Rules to Values

Applies a validation configuration's range rules to a numeric vector,
flagging out-of-range values.

## Usage

``` r
apply_validation_rules(values, property_name, validation_config)
```

## Arguments

- values:

  Numeric vector of values to validate

- property_name:

  Name of the property `values` represents

- validation_config:

  Validation configuration from
  [`create_validation_config()`](https://jjmaynard.github.io/soilSIM/reference/create_validation_config.md)

## Value

List with `violations` (logical vector, one per value) and
`rules_applied` (count of range rules evaluated)
