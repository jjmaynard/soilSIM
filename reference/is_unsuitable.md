# Identify Unsuitable Horizons

Identifies horizons that should be excluded from soil property modeling
based on horizon designation and characteristics. This is the
centralized function used by all modules.

## Usage

``` r
is_unsuitable(
  data,
  hzname_col = "hzname",
  strict_mode = TRUE,
  custom_exclusions = NULL,
  verbose = getOption("ssurgo.verbose", FALSE)
)
```

## Arguments

- data:

  Soil horizon data with horizon designation information

- hzname_col:

  Column name containing horizon designations (default = "hzname")

- strict_mode:

  Whether to use strict unsuitable criteria (default = TRUE)

- custom_exclusions:

  Additional horizon patterns to exclude

- verbose:

  Logical; if `TRUE`, temporarily raises the package's log level so
  `INFO`-level progress messages print for the duration of this call
  (default `FALSE` - quiet). See
  [`set_verbose_logging()`](https://jjmaynard.github.io/soilSIM/reference/set_verbose_logging.md).

## Value

Logical vector indicating unsuitable horizons
