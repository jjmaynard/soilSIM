# Standardize Property Names

Standardizes soil property names across different data sources using a
comprehensive mapping system.

## Usage

``` r
standardize_property_names(
  data,
  property_mapping = NULL,
  target_standard = "ssurgo",
  verbose = getOption("ssurgo.verbose", FALSE)
)
```

## Arguments

- data:

  Input data with potentially inconsistent property names

- property_mapping:

  Custom property mapping (optional)

- target_standard:

  Target naming standard ("ssurgo", "nrcs", "custom")

- verbose:

  Logical; if `TRUE`, temporarily raises the package's log level so
  `INFO`-level progress messages print for the duration of this call
  (default `FALSE` - quiet). See
  [`set_verbose_logging()`](https://jjmaynard.github.io/soilSIM/reference/set_verbose_logging.md).

## Value

Data with standardized property names
