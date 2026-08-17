# Setup Logging

Initializes comprehensive logging system for the workflow.

## Usage

``` r
setup_logging(
  log_file = NULL,
  log_level = "INFO",
  include_timestamp = TRUE,
  max_log_size = 10
)
```

## Arguments

- log_file:

  Log file path (NULL for console only)

- log_level:

  Log level ("DEBUG", "INFO", "WARN", "ERROR")

- include_timestamp:

  Whether to include timestamps

- max_log_size:

  Maximum log file size in MB

## Value

Logging configuration
