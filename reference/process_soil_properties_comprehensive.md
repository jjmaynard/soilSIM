# Comprehensive Soil Property Processing

The single property-infilling orchestrator. Infills a whole property set
in three phases - foundation (texture, bulk density, organic matter,
rock fragments), water retention, then the remaining chemical
properties - each property going through
[`infill_soil_property()`](https://jjmaynard.github.io/soilSIM/reference/infill_soil_property.md)'s
six-strategy hierarchy, with unsuitable horizons excluded throughout.
[`infill_soil_data()`](https://jjmaynard.github.io/soilSIM/reference/infill_soil_data.md)
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
