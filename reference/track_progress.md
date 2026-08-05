# Track Progress

Progress tracking utilities for long-running operations.

## Usage

``` r
track_progress(
  current,
  total,
  message = "Processing",
  update_frequency = 10,
  verbose = getOption("ssurgo.verbose", FALSE)
)
```

## Arguments

- current:

  Current progress value

- total:

  Total expected value

- message:

  Progress message

- update_frequency:

  Update frequency (every nth call)

- verbose:

  Logical; if `TRUE`, temporarily raises the package's log level so
  `INFO`-level progress messages print for the duration of this call
  (default `FALSE` - quiet). See
  [`set_verbose_logging()`](https://jjmaynard.github.io/soilSIM/reference/set_verbose_logging.md).
