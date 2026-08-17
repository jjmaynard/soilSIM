# Look up a compositional group's member property ids, in configured order

Reuses soilSIM's own `config$monte_carlo$composition_groups` convention
(see `R/distributions.R`'s
[`resolve_composition_groups()`](https://jjmaynard.github.io/soilSIM/reference/resolve_composition_groups.md))
rather than the original source bundle's separate
`PROPERTIES`/`config.R` global registry - order is trusted from the
config, not hardcoded to literal "clay"/"sand"/"silt" names, exactly
like the already-ported
[`fuse_texture_group_from_triplets()`](https://jjmaynard.github.io/soilSIM/reference/fuse_texture_group_from_triplets.md)
(`bayesian-updating.R`) already documents ("positional-role
placeholders... not an identity requirement").

## Usage

``` r
group_members(group, composition_groups)
```

## Arguments

- group:

  Composition group name (e.g. `"texture"`).

- composition_groups:

  A `config$monte_carlo$composition_groups`-shaped list, e.g.
  `list(texture = list(members = c("claytotal", "sandtotal", "silttotal")))`.

## Value

Character vector of member property ids, in configured order.
