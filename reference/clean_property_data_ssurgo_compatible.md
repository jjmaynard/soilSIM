# Clean Property Data (SSURGO Compatible)

Enhanced version incorporating SSURGO-specific cleaning logic while
using Module 8 utilities

## Usage

``` r
clean_property_data_ssurgo_compatible(
  df,
  property_name,
  validation_config = NULL,
  generate_report = TRUE,
  verbose = FALSE
)
```

## Arguments

- df:

  Input data frame

- property_name:

  Name of the property to clean

- validation_config:

  Optional validation configuration

- generate_report:

  Whether to generate detailed quality report

- verbose:

  Logical; provide progress messages

## Value

List containing cleaned data and optional quality report

## See also

[`advanced_string_parser_vectorized()`](https://jjmaynard.github.io/soilSIM/reference/advanced_string_parser_vectorized.md),
[`vectorized_type_conversion()`](https://jjmaynard.github.io/soilSIM/reference/vectorized_type_conversion.md),
and
[`apply_basic_range_limits()`](https://jjmaynard.github.io/soilSIM/reference/apply_basic_range_limits.md) -
the string-parsing/range-limit helpers this function calls are defined
in `data-infilling.R` (migrated from mod03), not in this file, after
consolidating this file's near-duplicate `*_working()`/`*_ssurgo()`
versions into those canonical implementations.
