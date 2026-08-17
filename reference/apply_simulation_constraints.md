# Apply Simulation Constraints (Enhanced)

Enhanced constraint application

## Usage

``` r
apply_simulation_constraints(
  simulation_results,
  properties,
  config,
  composition_plan = NULL,
  verbose = getOption("ssurgo.verbose", FALSE)
)
```

## Arguments

- simulation_results:

  Array of simulation results

- properties:

  Property names

- config:

  Configuration settings

- composition_plan:

  Optional result of
  [`resolve_composition_groups()`](https://jjmaynard.github.io/soilSIM/reference/resolve_composition_groups.md);
  when a group is active, its sum-to-100 constraint is already exact by
  construction (via
  [`restore_composition_properties()`](https://jjmaynard.github.io/soilSIM/reference/restore_composition_properties.md)'s
  ILR inverse), so
  [`get_sum_constraints()`](https://jjmaynard.github.io/soilSIM/reference/get_sum_constraints.md)
  skips it rather than re-rescaling.

- verbose:

  Logical; if `TRUE`, temporarily raises the package's log level so
  `INFO`-level progress messages print for the duration of this call
  (default `FALSE` - quiet). See
  [`set_verbose_logging()`](https://jjmaynard.github.io/soilSIM/reference/set_verbose_logging.md).

## Value

Constrained simulation results with quality metrics
