# Extract Simulation Parameters for One Horizon/Property (Enhanced)

Resolves the property's distribution family - from
`config$monte_carlo$property_distributions[[property_name]]$family`,
else `config$monte_carlo$distribution_type` (default `"triangular"`) -
and fits family-appropriate parameters from its SSURGO `_l/_r/_h`
triplet via `distributions.R`'s
[`fit_percentile_triplet()`](https://jjmaynard.github.io/soilSIM/reference/fit_percentile_triplet.md).
Previously this only ever returned `list(min=, mode=, max=)` regardless
of the requested `distribution_type`, so `"normal"`/`"beta"` silently
produced `NA` values downstream (their shape parameters were never
actually computed).

## Usage

``` r
extract_property_parameters(
  horizon_data,
  property_name,
  config,
  verbose = getOption("ssurgo.verbose", FALSE)
)
```

## Arguments

- horizon_data:

  One-row data frame/list for a single horizon.

- property_name:

  Property name (without `_l`/`_r`/`_h` suffix).

- config:

  Simulation configuration.

- verbose:

  Logical; if `TRUE`, temporarily raises the package's log level so
  `INFO`-level progress messages print for the duration of this call
  (default `FALSE` - quiet). See
  [`set_verbose_logging()`](https://jjmaynard.github.io/soilSIM/reference/set_verbose_logging.md).

## Value

`list(valid=, parameters=list(family=, fit=, source=))`.
