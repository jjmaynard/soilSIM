# Create a Validation Configuration

Creates an empty, pluggable validation-rule configuration that
[`add_range_rule()`](https://jjmaynard.github.io/soilSIM/reference/add_range_rule.md)
and
[`add_relationship_rule()`](https://jjmaynard.github.io/soilSIM/reference/add_relationship_rule.md)
can add rules to, and
[`apply_validation_rules()`](https://jjmaynard.github.io/soilSIM/reference/apply_validation_rules.md)
can apply.

## Usage

``` r
create_validation_config()
```

## Value

A validation configuration list with `range_rules`,
`relationship_rules`, `conditional_rules`, and `custom_rules` slots
