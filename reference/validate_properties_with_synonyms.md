# Validate Properties with Synonyms

Enhanced version that includes synonym matching

## Usage

``` r
validate_properties_with_synonyms(
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

  Property lookup source

- strict_mode:

  Enable strict validation

- max_invalid_pct:

  Maximum percentage of invalid properties allowed

- performance_threshold:

  Performance warning threshold

## Value

Enhanced validation results with synonym matching
