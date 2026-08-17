# Assess Group Representation

Assesses whether groups (`cokey`) in `simulation_data` are adequately
represented (each has at least `criteria$min_samples_per_group` rows).

## Usage

``` r
assess_group_representation(simulation_data, criteria)
```

## Arguments

- simulation_data:

  Simulation data with a `cokey` column.

- criteria:

  List, optionally with `min_samples_per_group` (default 10).

## Value

List with `groups_well_represented` (fraction),
`underrepresented_groups`, `representation_adequate`.
