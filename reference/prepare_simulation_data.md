# Prepare Simulation Data (Enhanced)

Enhanced data preparation

## Usage

``` r
prepare_simulation_data(
  soil_data,
  properties,
  config,
  composition_plan = NULL,
  verbose = getOption("ssurgo.verbose", FALSE)
)
```

## Arguments

- soil_data:

  Input soil data

- properties:

  Properties to prepare (may include composition-group
  pseudo-properties, e.g. "ilr1"/"ilr2" - see `composition_plan`)

- config:

  Configuration settings

- composition_plan:

  Optional result of
  [`resolve_composition_groups()`](https://jjmaynard.github.io/soilSIM/reference/resolve_composition_groups.md);
  when supplied, data-availability checks are resolved back to the real
  underlying SSURGO columns (e.g.
  `sandtotal_r`/`claytotal_r`/`silttotal_r`) for any pseudo-property in
  `properties`, since a literal `"ilr1_r"` column never exists.

- verbose:

  Logical; if `TRUE`, temporarily raises the package's log level so
  `INFO`-level progress messages print for the duration of this call
  (default `FALSE` - quiet). See
  [`set_verbose_logging()`](https://jjmaynard.github.io/soilSIM/reference/set_verbose_logging.md).

## Value

Prepared simulation data
