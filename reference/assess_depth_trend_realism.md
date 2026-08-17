# Assess Depth Trend Realism

Enhanced depth trend validation using Module 5 prediction functions and
Module 0 utilities.

## Usage

``` r
assess_depth_trend_realism(
  gp_models,
  criteria = NULL,
  verbose = getOption("ssurgo.verbose", FALSE)
)
```

## Arguments

- gp_models:

  GP models

- criteria:

  Realism criteria

- verbose:

  Logical; if `TRUE`, temporarily raises the package's log level so
  `INFO`-level progress messages print for the duration of this call
  (default `FALSE` - quiet). See
  [`set_verbose_logging()`](https://jjmaynard.github.io/soilSIM/reference/set_verbose_logging.md).

## Value

Depth trend realism assessment
