# Read Soil Data

Robust soil data reading with automatic format detection and validation.

## Usage

``` r
read_soil_data(
  file_path,
  file_type = "auto",
  validate_on_load = TRUE,
  sheet_name = NULL,
  verbose = getOption("ssurgo.verbose", FALSE)
)
```

## Arguments

- file_path:

  Path to data file

- file_type:

  File type ("auto", "csv", "xlsx", "rds", "txt")

- validate_on_load:

  Whether to validate data after loading

- sheet_name:

  Sheet name for Excel files

- verbose:

  Logical; if `TRUE`, temporarily raises the package's log level so
  `INFO`-level progress messages print for the duration of this call
  (default `FALSE` - quiet). See
  [`set_verbose_logging()`](https://jjmaynard.github.io/soilSIM/reference/set_verbose_logging.md).

## Value

Loaded soil data

## Usage note

Fully implemented and exported, but not currently called from anywhere
else in the package outside a `@seealso` doc link (confirmed by source
grep) - standalone public API for callers who need it, not dead/broken
code. See also
[`write_soil_data`](https://jjmaynard.github.io/soilSIM/reference/write_soil_data.md).
