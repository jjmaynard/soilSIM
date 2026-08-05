# Write Soil Data

Standardized soil data writing with metadata and backup options.

## Usage

``` r
write_soil_data(
  data,
  file_path,
  file_type = "csv",
  include_metadata = TRUE,
  create_backup = TRUE,
  verbose = getOption("ssurgo.verbose", FALSE)
)
```

## Arguments

- data:

  Data to write

- file_path:

  Output file path

- file_type:

  Output file type ("csv", "xlsx", "rds")

- include_metadata:

  Whether to include metadata

- create_backup:

  Whether to create backup

- verbose:

  Logical; if `TRUE`, temporarily raises the package's log level so
  `INFO`-level progress messages print for the duration of this call
  (default `FALSE` - quiet). See
  [`set_verbose_logging()`](https://jjmaynard.github.io/soilSIM/reference/set_verbose_logging.md).

## Value

Success status
