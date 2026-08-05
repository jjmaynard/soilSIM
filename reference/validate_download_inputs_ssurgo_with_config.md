# Enhanced SSURGO Download Input Validation (Wrapper)

Convenience wrapper that provides additional configuration options and
integrates with Module 8 configuration management

## Usage

``` r
validate_download_inputs_ssurgo_with_config(
  aoi_wkt,
  properties,
  include_restrictions = TRUE,
  config_file = NULL,
  log_validation = TRUE
)
```

## Arguments

- aoi_wkt:

  Area of interest in WKT format

- properties:

  Vector of property names

- include_restrictions:

  Logical for restriction data inclusion

- config_file:

  Optional configuration file path

- log_validation:

  Whether to enable detailed logging

## Value

Enhanced validation results
