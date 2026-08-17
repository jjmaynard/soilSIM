# Validate Soil Properties (Simple Generic Version)

Simple, extensible property validation that works with any property
source. Focuses on core validation logic without excessive complexity.

## Usage

``` r
validate_properties(
  properties,
  property_lookup = "ssurgo",
  strict_mode = TRUE,
  max_invalid_pct = 50,
  performance_threshold = 10
)
```

## Arguments

- properties:

  Vector of property names to validate

- property_lookup:

  Property lookup source (function, vector, data.frame, or source name)

- strict_mode:

  Enable strict validation (default: TRUE)

- max_invalid_pct:

  Maximum percentage of invalid properties allowed (default: 50)

- performance_threshold:

  Warn if more than this many properties (default: 10)

## Value

Simple validation results list
