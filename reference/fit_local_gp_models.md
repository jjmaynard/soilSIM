# Fit Local GP Models

Enhanced version with Module 8 validation and configuration management.

## Usage

``` r
fit_local_gp_models(
  cokey_data,
  properties,
  config = NULL,
  gp_control = c(20, 10, 2),
  verbose = getOption("ssurgo.verbose", FALSE)
)
```

## Arguments

- cokey_data:

  Simulation data for a single cokey

- properties:

  Properties to model

- config:

  Configuration settings

- gp_control:

  Passed through to
  [`fit_local_gp_model_single()`](https://jjmaynard.github.io/soilSIM/reference/fit_local_gp_model_single.md)'s
  `gp_control` - see its docs for why the default is much smaller than
  [`GPfit::GP_fit()`](https://rdrr.io/pkg/GPfit/man/GP_fit.html)'s own
  default.

- verbose:

  Logical; if `TRUE`, temporarily raises the package's log level so
  `INFO`-level progress messages print for the duration of this call
  (default `FALSE` - quiet). See
  [`set_verbose_logging()`](https://jjmaynard.github.io/soilSIM/reference/set_verbose_logging.md).

## Value

List of fitted local GP models
