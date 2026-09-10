# Add a Relationship Validation Rule

Adds a cross-property relationship rule (e.g. a sum-to-100 texture
constraint) to a validation configuration. Relationship rules are
recorded on the config but are not consumed by
[`apply_validation_rules()`](https://jjmaynard.github.io/soilSIM/reference/apply_validation_rules.md),
which enforces only `range_rules`.

## Usage

``` r
add_relationship_rule(
  config,
  properties,
  relationship_type,
  expected_sum = NULL,
  tolerance = 0.1
)
```

## Arguments

- config:

  Validation configuration from
  [`create_validation_config()`](https://jjmaynard.github.io/soilSIM/reference/create_validation_config.md)

- properties:

  Character vector of properties involved in the relationship

- relationship_type:

  Type of relationship, e.g. `"sum"`

- expected_sum:

  Expected sum, when `relationship_type == "sum"`

- tolerance:

  Allowed tolerance around `expected_sum`

## Value

Updated validation configuration
