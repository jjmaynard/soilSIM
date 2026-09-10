# Validate a Horizon's Fitted Distribution Parameters

Thin wrapper around `distributions.R`'s
[`validate_fit_parameters()`](https://jjmaynard.github.io/soilSIM/reference/validate_fit_parameters.md).

## Usage

``` r
validate_distribution_parameters(parameters, distribution_type, property_name)
```

## Arguments

- parameters:

  A params object from
  [`extract_property_parameters()`](https://jjmaynard.github.io/soilSIM/reference/extract_property_parameters.md)
  (`list(family=, fit=, source=)`).

- distribution_type:

  Fallback family if `parameters$family` is absent.

- property_name:

  Unused; kept for interface compatibility with callers.

## Value

`list(valid=, message=)`.
