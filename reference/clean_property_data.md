# Property Data Cleaning

The single property-value cleaner for the whole package: parse
string/factor cells to numeric, null non-finite values, optionally flag
statistical outliers, then clamp to hardcoded physical-plausibility
bounds. Both the tabular pipeline
([`process_ssurgo_data()`](https://jjmaynard.github.io/soilSIM/reference/process_ssurgo_data.md))
and the infilling hierarchy
([`infill_soil_property()`](https://jjmaynard.github.io/soilSIM/reference/infill_soil_property.md))
call this.

## Usage

``` r
clean_property_data(
  df,
  property_name,
  outlier_policy = c("soil_aware", "aggressive_iqr", "none"),
  validation_config = NULL,
  generate_report = FALSE,
  verbose = FALSE
)
```

## Arguments

- df:

  Input data frame.

- property_name:

  Name of the property to clean.

- outlier_policy:

  How to flag statistical outliers on the `_r` column, before the
  physical-bound clamp:

  `"soil_aware"` (default)

  :   [`detect_statistical_outliers_soil_aware()`](https://jjmaynard.github.io/soilSIM/reference/detect_statistical_outliers_soil_aware.md) -
      flags only physically impossible values for pH / texture / bulk
      density, and a very conservative `IQR x 5` for everything else.
      The right default for a package whose job is uncertainty
      propagation: a real tail value is signal, not noise.

  `"aggressive_iqr"`

  :   generic `IQR x 3` for every property alike. Opt-in only.

  `"none"`

  :   skip outlier flagging; still parse, null non-finite, and clamp.

- validation_config:

  Optional per-property range rule (a `list(min=, max=)`-shaped entry);
  violations are counted into the report, values are not changed.

- generate_report:

  Whether to return a cleaning report alongside the data.

- verbose:

  Whether to provide progress messages.

## Value

`list(data = <cleaned df>)`, plus `report` when
`generate_report = TRUE`.
