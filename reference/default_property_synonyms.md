# Get Default Synonyms for Property Sources

Returns common synonym mappings for different property sources

## Usage

``` r
default_property_synonyms(property_lookup)
```

## Arguments

- property_lookup:

  Property lookup source

## Value

Named list of synonyms

## Usage note

Fully implemented and exported, but not currently called from anywhere
else in the package (confirmed by source grep) - standalone public API
for callers who need it, not dead/broken code. Its `"ssurgo"` branch
shares the same canonical synonym data as
`get_default_property_mapping()` (internal) and
[`validate_properties_with_synonyms`](https://jjmaynard.github.io/soilSIM/reference/validate_properties_with_synonyms.md)
(see `.ssurgo_property_synonyms()`).

## See also

Other utilities:
[`available_properties()`](https://jjmaynard.github.io/soilSIM/reference/available_properties.md),
[`default_config()`](https://jjmaynard.github.io/soilSIM/reference/default_config.md),
[`default_diagnostics_config()`](https://jjmaynard.github.io/soilSIM/reference/default_diagnostics_config.md),
[`is_unsuitable()`](https://jjmaynard.github.io/soilSIM/reference/is_unsuitable.md),
[`predefined_properties()`](https://jjmaynard.github.io/soilSIM/reference/predefined_properties.md),
[`validate_data_quality()`](https://jjmaynard.github.io/soilSIM/reference/validate_data_quality.md)
