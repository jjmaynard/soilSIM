# Prepare Simulation Parameters (Enhanced)

Enhanced parameter preparation. Composition-group pseudo-properties
(e.g. `"ilr1"`/`"ilr2"`) are special-cased: BOTH are fit together, once
per horizon, via `distributions.R`'s
[`estimate_ilr_moments_mc()`](https://jjmaynard.github.io/soilSIM/reference/estimate_ilr_moments_mc.md)
from the group's real member `_l/_r/_h` triplets - their joint
covariance only makes sense computed jointly, not by looping
[`extract_property_parameters()`](https://jjmaynard.github.io/soilSIM/reference/extract_property_parameters.md)
per pseudo name. The Monte Carlo covariance's ilr1\<-\>ilr2 OFF-DIAGONAL
is intentionally discarded (only the marginal SDs are kept): the
cross-horizon empirical correlation-matrix step
([`estimate_property_correlations()`](https://jjmaynard.github.io/soilSIM/reference/estimate_property_correlations.md))
owns ilr1\<-\>ilr2 correlation instead, so the two sources of
correlation (per-observation uncertainty propagation vs. cross-sample
empirical correlation) are never double-counted.

## Usage

``` r
prepare_simulation_parameters(
  simulation_data,
  properties,
  config,
  composition_plan = NULL,
  verbose = getOption("ssurgo.verbose", FALSE)
)
```

## Arguments

- simulation_data:

  Filtered simulation data

- properties:

  Properties to extract parameters for (may include composition-group
  pseudo-properties)

- config:

  Configuration settings

- composition_plan:

  Optional result of
  [`resolve_composition_groups()`](https://jjmaynard.github.io/soilSIM/reference/resolve_composition_groups.md).

- verbose:

  Logical; if `TRUE`, temporarily raises the package's log level so
  `INFO`-level progress messages print for the duration of this call
  (default `FALSE` - quiet). See
  [`set_verbose_logging()`](https://jjmaynard.github.io/soilSIM/reference/set_verbose_logging.md).

## Value

Enhanced simulation parameters
