# Validate Download Inputs (SSURGO-specific)

SSURGO-specific input validation that uses Module 8 general validation
functions for parameter validation, logging, and error handling, with
custom SSURGO-specific business logic for WKT geometry and property
validation.

## Usage

``` r
validate_download_inputs_ssurgo(
  aoi_wkt,
  properties,
  include_restrictions = TRUE,
  strict_geometry = TRUE,
  max_area_deg2 = 100
)
```

## Arguments

- aoi_wkt:

  Area of interest in WKT format

- properties:

  Vector of property names

- include_restrictions:

  Logical for restriction data inclusion

- strict_geometry:

  Whether to use strict geometry validation

- max_area_deg2:

  Maximum allowed area in square degrees

## Value

List with validation results including Module 8 metadata
