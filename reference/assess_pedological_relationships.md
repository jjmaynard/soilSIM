# Assess Pedological Relationships

Checks a small set of cheap, well-established pedological relationship
directions when the relevant properties are present: clay content
positively associated with CEC, and organic matter negatively associated
with depth.

## Usage

``` r
assess_pedological_relationships(simulation_data, criteria)
```

## Arguments

- simulation_data:

  Simulation data.

- criteria:

  Unused (no per-relationship configuration currently defined upstream).

## Value

List with `texture_relationship_quality` (fraction of rows with
sand+clay+silt within 5 of 100), `chemical_relationship_quality`
(rescaled clay-CEC correlation), `physical_relationship_quality`
(rescaled OM-depth correlation) - each `NA` when the needed properties
are absent.
