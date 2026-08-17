# Collapse simulated ILR pseudo-properties back to raw composition members

For each active composition group, applies
[`ilr_inverse()`](https://jjmaynard.github.io/soilSIM/reference/ilr_inverse.md)
per horizon (vectorized over realizations) and reassembles a
`[horizon, property, realization]` array dimnamed over the ORIGINAL
`properties` vector - non-group properties are copied through unchanged,
group members are filled from the inverse-ILR output (guaranteed
sum-to-`total` and in-bounds for every realization by construction), and
the pseudo-properties are dropped.

## Usage

``` r
restore_composition_properties(
  simulation_results,
  sim_properties,
  properties,
  groups
)
```

## Arguments

- simulation_results:

  A `[horizon, property, realization]` array dimnamed over
  `sim_properties`.

- sim_properties:

  The (possibly ILR-substituted) property vector `simulation_results` is
  dimnamed over.

- properties:

  The caller-facing, original property vector.

- groups:

  The `groups` element of
  [`resolve_composition_groups()`](https://jjmaynard.github.io/soilSIM/reference/resolve_composition_groups.md)'s
  return value.

## Value

An array dimnamed over `properties`.
