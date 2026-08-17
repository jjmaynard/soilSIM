# Resolve composition-group pseudo-properties to their real backing columns

Resolve composition-group pseudo-properties to their real backing
columns

## Usage

``` r
resolve_real_properties(properties, composition_plan = NULL)
```

## Arguments

- properties:

  Character vector, possibly including pseudo-properties (e.g.
  "ilr1"/"ilr2") from an active composition group.

- composition_plan:

  Optional result of
  [`resolve_composition_groups()`](https://jjmaynard.github.io/soilSIM/reference/resolve_composition_groups.md).

## Value

Character vector of real property names (pseudo-properties expanded to
their group's `members`; everything else passed through unchanged).
