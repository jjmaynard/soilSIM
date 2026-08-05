# Determine sum-to-100 constraints (e.g. texture) to apply post-simulation

Determine sum-to-100 constraints (e.g. texture) to apply post-simulation

## Usage

``` r
get_sum_constraints(properties, composition_plan = NULL)
```

## Arguments

- properties:

  Character vector of property names.

- composition_plan:

  Optional result of
  [`resolve_composition_groups()`](https://jjmaynard.github.io/soilSIM/reference/resolve_composition_groups.md);
  when the texture group is active, sum-to-100 is already exact by
  construction (via
  [`restore_composition_properties()`](https://jjmaynard.github.io/soilSIM/reference/restore_composition_properties.md)'s
  ILR inverse), so no proportional-rescale entry is returned for it -
  re-rescaling would distort the correlation structure the ILR path
  exists to preserve.
