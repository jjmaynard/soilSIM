# Comprehensive Soil Property Processing

Complete soil property processing workflow that handles multiple
properties with automatic exclusion of unsuitable horizons and
intelligent infilling strategies.

## Usage

``` r
process_soil_properties_comprehensive(
  df,
  properties = NULL,
  max_depth = 250,
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
