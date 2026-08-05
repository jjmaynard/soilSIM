# Add a Range Validation Rule

Adds a rule to a validation configuration flagging values of `property`
outside `[min_val, max_val]`.

## Usage

``` r
add_range_rule(config, property, min_val, max_val, severity = "error")
```

## Arguments

- config:

  Validation configuration from
  [`create_validation_config()`](https://jjmaynard.github.io/soilSIM/reference/create_validation_config.md)

- property:

  Name of the property this rule applies to

- min_val:

  Minimum acceptable value

- max_val:

  Maximum acceptable value

- severity:

  Rule severity, e.g. `"error"` or `"warning"`

## Value

Updated validation configuration
