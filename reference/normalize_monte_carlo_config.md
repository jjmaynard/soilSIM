# Normalize a Monte Carlo `simulation_config` to be Properly Nested

[`generate_monte_carlo_realizations()`](https://jjmaynard.github.io/soilSIM/reference/generate_monte_carlo_realizations.md)'s
own `@examples` have always shown a FLAT `simulation_config` (e.g.
`list(distribution_type = "normal", max_depth = 200)`), but
[`get_monte_carlo_defaults()`](https://jjmaynard.github.io/soilSIM/reference/get_monte_carlo_defaults.md)
nests every Monte Carlo setting under `$monte_carlo`, and
[`merge_configurations()`](https://jjmaynard.github.io/soilSIM/reference/merge_configurations.md)
merges strictly by matching key path - so a flat user config silently
landed as a new, unused top-level `config$distribution_type` and had NO
effect at all (confirmed by testing: `distribution_type = "normal"`
passed flat left `config$monte_carlo$distribution_type` at its default
`"triangular"`, every time). This accepts either shape: already-nested
(`list(monte_carlo = list(...))`) configs pass through unchanged; a flat
config whose names match known `monte_carlo` settings gets wrapped.

## Usage

``` r
normalize_monte_carlo_config(
  simulation_config,
  default_config,
  verbose = getOption("ssurgo.verbose", FALSE)
)
```

## Arguments

- simulation_config:

  User-supplied config (flat or nested).

- default_config:

  Result of
  [`get_monte_carlo_defaults()`](https://jjmaynard.github.io/soilSIM/reference/get_monte_carlo_defaults.md),
  used only to recognize which flat key names belong under
  `monte_carlo`.

- verbose:

  Logical; if `TRUE`, temporarily raises the package's log level so
  `INFO`-level progress messages print for the duration of this call
  (default `FALSE` - quiet). See
  [`set_verbose_logging()`](https://jjmaynard.github.io/soilSIM/reference/set_verbose_logging.md).

## Value

`simulation_config`, wrapped under `monte_carlo` if it was flat.
