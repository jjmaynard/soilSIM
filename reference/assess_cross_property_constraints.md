# Assess Cross-Property Constraints

Checks the texture sum-to-100 constraint (when `sandtotal`/`claytotal`/
`silttotal` are present) and any explicit `cross_property_rules` (in the
same shape as
[`create_validation_config()`](https://jjmaynard.github.io/soilSIM/reference/create_validation_config.md)'s
`relationship_rules` -
`list(properties=, type=, expected_sum=, tolerance=)`).

## Usage

``` r
assess_cross_property_constraints(simulation_data, cross_property_rules)
```

## Arguments

- simulation_data:

  Simulation data.

- cross_property_rules:

  List of relationship rules, or `NULL`.

## Value

List with `texture_sum_violations`, `relationship_violations`,
`overall_compliance`.
