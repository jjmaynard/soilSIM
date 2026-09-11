# Clean Property Data (SSURGO Compatible) - deprecated

**Deprecated.** Use
[`clean_property_data()`](https://jjmaynard.github.io/soilSIM/reference/clean_property_data.md)
directly. This shim forwards to it with
`outlier_policy = "aggressive_iqr"` (a generic `IQR x 3`) and
`generate_report = TRUE`. New code should choose an `outlier_policy`
explicitly; `"soil_aware"` is the recommended default.

## Usage

``` r
clean_ssurgo_property_data(
  df,
  property_name,
  validation_config = NULL,
  generate_report = TRUE,
  verbose = FALSE
)
```

## Arguments

- df:

  Input data frame.

- property_name:

  Name of the property to clean.

- validation_config:

  Optional per-property range rule (a `list(min=, max=)`-shaped entry);
  violations are counted into the report, values are not changed.

- generate_report:

  Whether to return a cleaning report alongside the data.

- verbose:

  Whether to provide progress messages.

## Value

List containing cleaned data and (by default) a quality report.

## See also

[`clean_property_data()`](https://jjmaynard.github.io/soilSIM/reference/clean_property_data.md)
