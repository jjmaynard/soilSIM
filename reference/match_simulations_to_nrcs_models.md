# Match Simulations to NRCS Models

Looks up the NRCS GP model group for a cokey in cokey_mapping, returning
fallback_group when no mapping is found.

## Usage

``` r
match_simulations_to_nrcs_models(
  cokey,
  cokey_mapping,
  fallback_group = "general_pool",
  verbose = getOption("ssurgo.verbose", FALSE)
)
```

## Arguments

- cokey:

  Target cokey

- cokey_mapping:

  Mapping from match_soils_to_gp_models()

- fallback_group:

  Default group if matching fails

- verbose:

  Logical; if `TRUE`, temporarily raises the package's log level so
  `INFO`-level progress messages print for the duration of this call
  (default `FALSE` - quiet). See
  [`set_verbose_logging()`](https://jjmaynard.github.io/soilSIM/reference/set_verbose_logging.md).

## Value

Model group name
