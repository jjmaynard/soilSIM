# Get Validation Defaults by Context

Get appropriate default validation parameters for different contexts

## Usage

``` r
default_diagnostics_config(context, crs)
```

## Arguments

- context:

  Validation context

- crs:

  Coordinate reference system

## Value

Default validation parameters

## Usage note

Fully implemented and exported, but not currently called from anywhere
else in the package (confirmed by source grep) - standalone public API,
not dead/broken code.
[`validate_wkt_geometry`](https://jjmaynard.github.io/soilSIM/reference/validate_wkt_geometry.md)
(the package's live WKT validation entry point) hardcodes its own
context-appropriate defaults directly in its function signature rather
than calling this at runtime, so the two default sets should be kept in
sync by hand if either changes.

## See also

Other utilities:
[`available_properties()`](https://jjmaynard.github.io/soilSIM/reference/available_properties.md),
[`default_config()`](https://jjmaynard.github.io/soilSIM/reference/default_config.md),
[`default_property_synonyms()`](https://jjmaynard.github.io/soilSIM/reference/default_property_synonyms.md),
[`is_unsuitable()`](https://jjmaynard.github.io/soilSIM/reference/is_unsuitable.md),
[`predefined_properties()`](https://jjmaynard.github.io/soilSIM/reference/predefined_properties.md),
[`validate_data_quality()`](https://jjmaynard.github.io/soilSIM/reference/validate_data_quality.md)
