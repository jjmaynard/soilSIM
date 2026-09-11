# Comprehensive Soil Property Processing

The single property-infilling orchestrator. Infills a whole property set
in three phases - foundation (texture, bulk density, organic matter,
rock fragments), water retention, then the remaining chemical
properties - each property going through
[`infill_soil_property()`](https://jjmaynard.github.io/soilSIM/reference/infill_soil_property.md)'s
six-strategy hierarchy, with unsuitable horizons excluded throughout.
[`infill_ssurgo_data()`](https://jjmaynard.github.io/soilSIM/reference/infill_ssurgo_data.md)
is a thin wrapper that calls this with the standard SSURGO property set.

## Usage

``` r
process_soil_properties_comprehensive(
  df,
  properties = NULL,
  max_depth = DEFAULT_MAX_DEPTH_CM,
  water_retention_method = c("saxton_rawls", "generic"),
  remove_unsuitable = FALSE,
  remove_incomplete = FALSE,
  required_properties = NULL,
  verbose = getOption("ssurgo.verbose", FALSE)
)
```

## Arguments

- df:

  Input soil data frame

- properties:

  Vector of properties to process (NULL = auto-detect)

- max_depth:

  Maximum depth for processing (default: 250 cm)

- water_retention_method:

  How to fill `wthirdbar`/`wfifteenbar`: `"saxton_rawls"` (default) runs
  the Saxton-Rawls pedotransfer function where texture + bulk density
  allow and the generic hierarchy elsewhere; `"generic"` sends water
  retention through the six-strategy hierarchy like any other property
  (which itself falls back to Saxton-Rawls at strategy 5).

- remove_unsuitable:

  Whether to remove unsuitable horizons from output

- remove_incomplete:

  Whether to remove incomplete rows

- required_properties:

  Vector of properties that must be complete

- verbose:

  Whether to print detailed progress messages

## Value

Data frame with processed and infilled properties

## See also

Other ssurgo-adapter:
[`add_range_rule()`](https://jjmaynard.github.io/soilSIM/reference/add_range_rule.md),
[`add_relationship_rule()`](https://jjmaynard.github.io/soilSIM/reference/add_relationship_rule.md),
[`apply_validation_rules()`](https://jjmaynard.github.io/soilSIM/reference/apply_validation_rules.md),
[`clean_property_data()`](https://jjmaynard.github.io/soilSIM/reference/clean_property_data.md),
[`create_custom_property_config()`](https://jjmaynard.github.io/soilSIM/reference/create_custom_property_config.md),
[`create_validation_config()`](https://jjmaynard.github.io/soilSIM/reference/create_validation_config.md),
[`default_property_config()`](https://jjmaynard.github.io/soilSIM/reference/default_property_config.md),
[`impute_rfv_values()`](https://jjmaynard.github.io/soilSIM/reference/impute_rfv_values.md),
[`infill_missing_property_data()`](https://jjmaynard.github.io/soilSIM/reference/infill_missing_property_data.md),
[`infill_water_retention_saxton_rawls_integrated()`](https://jjmaynard.github.io/soilSIM/reference/infill_water_retention_saxton_rawls_integrated.md),
[`learn_property_ranges()`](https://jjmaynard.github.io/soilSIM/reference/learn_property_ranges.md),
[`process_ssurgo_components()`](https://jjmaynard.github.io/soilSIM/reference/process_ssurgo_components.md),
[`process_ssurgo_data()`](https://jjmaynard.github.io/soilSIM/reference/process_ssurgo_data.md),
[`process_ssurgo_horizons()`](https://jjmaynard.github.io/soilSIM/reference/process_ssurgo_horizons.md),
[`property_contextual_ranges()`](https://jjmaynard.github.io/soilSIM/reference/property_contextual_ranges.md),
[`summarize_unsuitable_horizons()`](https://jjmaynard.github.io/soilSIM/reference/summarize_unsuitable_horizons.md)
