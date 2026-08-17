# Estimate a Correlation Matrix from Simulation Parameters (Real Implementation)

Builds one representative-value-per-horizon data frame - the per-horizon
fitted mean/mode for each property (for composition-group
pseudo-properties like ilr1/ilr2 there's no raw `_r` column, so the
fitted mean IS the representative value) - plus `genhz` from
`simulation_data` if present, then calls `distributions.R`'s
[`estimate_correlation_matrix_robust()`](https://jjmaynard.github.io/soilSIM/reference/estimate_correlation_matrix_robust.md)
([`Hmisc::rcorr()`](https://rdrr.io/pkg/Hmisc/man/rcorr.html) +
PD-repair, genhz-stratified when possible). Previously this was a
placeholder always returning `valid=FALSE`, so `auto_correlation=TRUE`
silently fell back to the identity matrix.

## Usage

``` r
estimate_property_correlations(
  simulation_params,
  properties,
  config,
  simulation_data = NULL,
  verbose = getOption("ssurgo.verbose", FALSE)
)
```

## Arguments

- simulation_params:

  List of per-horizon parameter lists.

- properties:

  Character vector of property names.

- config:

  Simulation configuration.

- simulation_data:

  Optional raw per-horizon data frame (same row order as
  `simulation_params`); its `genhz` column, if present, enables
  genhz-stratified estimation (an explicit `genhz` column always takes
  precedence over auto-derivation from `hzname`).

- verbose:

  Logical; if `TRUE`, temporarily raises the package's log level so
  `INFO`-level progress messages print for the duration of this call
  (default `FALSE` - quiet). See
  [`set_verbose_logging()`](https://jjmaynard.github.io/soilSIM/reference/set_verbose_logging.md).

## Value

`list(valid=, correlation_matrix=, info=list(method=, n_obs=))`.

## Details

When `config$monte_carlo$correlation_fallback == "kssl_global"` (opt-in;
defaults to `"identity"`, today's exact behavior), a group that fails
empirical estimation falls back to that group's own KSSL reference
correlation matrix (`kssl-reference-correlations.R`) instead of being
dropped, and the final identity fallback becomes a KSSL-pooled matrix
instead. `genhz` is auto-derived from `simulation_data$hzname` via
[`classify_genhz()`](https://jjmaynard.github.io/soilSIM/reference/classify_genhz.md)
when this is requested and `simulation_data$genhz` isn't already
present - scoped narrowly behind the opt-in flag so a caller who has
`hzname` but never requests `"kssl_global"` sees no behavior change at
all.
