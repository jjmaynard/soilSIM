# Validate a Property Configuration

Validates that a property configuration list has the required fields and
that they are well-formed.

## Usage

``` r
validate_property_config(config, property_name)
```

## Arguments

- config:

  Property configuration list

- property_name:

  Name of the property (for error messages)

## Value

`TRUE` if valid; stops with an error if invalid
