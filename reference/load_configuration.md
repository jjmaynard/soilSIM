# Load Configuration

Loads configuration from JSON or R files with validation and defaults.

## Usage

``` r
load_configuration(
  config_path,
  config_type = "auto",
  validate_config = TRUE,
  verbose = getOption("ssurgo.verbose", FALSE)
)
```

## Arguments

- config_path:

  Path to configuration file

- config_type:

  Configuration type ("json", "r", "auto")

- validate_config:

  Whether to validate configuration

- verbose:

  Logical; if `TRUE`, temporarily raises the package's log level so
  `INFO`-level progress messages print for the duration of this call
  (default `FALSE` - quiet). See
  [`set_verbose_logging()`](https://jjmaynard.github.io/soilSIM/reference/set_verbose_logging.md).

## Value

Loaded configuration
