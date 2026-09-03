# Get Validation Defaults by Context

Get appropriate default validation parameters for different contexts

## Usage

``` r
get_validation_defaults(context, crs)
```

## Arguments

- context:

  Validation context

- crs:

  Coordinate reference system

## Value

Default validation parameters

## Usage note

Fully implemented and exported, but not currently called from anywhere
else in the package (confirmed by source grep) - standalone public API,
not dead/broken code.
[`validate_wkt_geometry`](https://jjmaynard.github.io/soilSIM/reference/validate_wkt_geometry.md)
(the package's live WKT validation entry point) hardcodes its own
context-appropriate defaults directly in its function signature rather
than calling this at runtime, so the two default sets should be kept in
sync by hand if either changes.
