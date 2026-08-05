# Handle Missing Values

Comprehensive missing value handling with multiple strategies.

## Usage

``` r
handle_missing_values(
  data,
  strategy = "interpolate",
  columns = NULL,
  group_by = NULL,
  verbose = getOption("ssurgo.verbose", FALSE)
)
```

## Arguments

- data:

  Input data

- strategy:

  Missing value strategy ("remove", "interpolate", "mean", "median",
  "forward_fill")

- columns:

  Columns to process (NULL = all numeric columns)

- group_by:

  Grouping columns for grouped imputation

- verbose:

  Logical; if `TRUE`, temporarily raises the package's log level so
  `INFO`-level progress messages print for the duration of this call
  (default `FALSE` - quiet). See
  [`set_verbose_logging()`](https://jjmaynard.github.io/soilSIM/reference/set_verbose_logging.md).

## Value

Data with handled missing values
