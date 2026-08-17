# Expand a requested properties vector around active composition groups

If ALL members of a configured composition group (e.g. `"texture"`:
clay/sand/silt) are present in `properties`, they are replaced (position
of the first member preserved) with the group's pseudo-property names
(e.g. `"ilr1"`/`"ilr2"`) so they can flow through the existing generic
Cholesky-copula simulation machinery like any other property. If only
some members are present, the group stays inactive (a logged WARN) and
those properties simulate independently via the legacy path.

## Usage

``` r
resolve_composition_groups(properties, config)
```

## Arguments

- properties:

  Character vector of requested property names.

- config:

  A Monte Carlo config list; reads
  `config$monte_carlo$composition_groups` (a named list, each entry
  `list(members=, pseudo=)`).

## Value

`list(sim_properties=, groups=)` where `groups` is a named list of
`list(members=, pseudo=, active=)`.

## Details

`members` declares which real property occupies
[`ilr_forward()`](https://jjmaynard.github.io/soilSIM/reference/ilr_forward.md)/
[`ilr_inverse()`](https://jjmaynard.github.io/soilSIM/reference/ilr_inverse.md)'s
position 1/2/3 (their `clay`/`sand`/`silt` naming is a positional-role
placeholder, not an identity requirement - see the header comment in the
ILR section of this file).
[`restore_composition_properties()`](https://jjmaynard.github.io/soilSIM/reference/restore_composition_properties.md)
maps
[`ilr_inverse()`](https://jjmaynard.github.io/soilSIM/reference/ilr_inverse.md)'s
fixed `clay`/`sand`/`silt` output columns onto
`members[1]`/`members[2]`/`members[3]` by POSITION, not by name
matching.
