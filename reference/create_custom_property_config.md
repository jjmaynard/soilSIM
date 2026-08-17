# Create Custom Property Configuration

Creates a custom configuration for specialized soil properties.

## Usage

``` r
create_custom_property_config(
  property_name,
  property_type = "generic",
  units = "unknown",
  typical_range = NULL,
  fallback_range = 5,
  related_properties = NULL,
  special_options = NULL
)
```

## Arguments

- property_name:

  Name of the property

- property_type:

  Type of property

- units:

  Units of measurement

- typical_range:

  Numeric vector c(min, max) of typical values

- fallback_range:

  Default spread for range calculation

- related_properties:

  Character vector of related property names

- special_options:

  List of additional options

## Value

Property configuration list
