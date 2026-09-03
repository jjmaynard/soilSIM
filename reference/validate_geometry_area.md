# Validate Geometry Area

Validate geometry area against specified limits

## Usage

``` r
validate_geometry_area(geom, area_limits, context = "general")
```

## Arguments

- geom:

  Terra geometry object

- area_limits:

  List with min/max area limits

- context:

  Validation context for appropriate messaging

## Value

Area validation results

## Usage note

Fully implemented and exported, but not currently called from anywhere
else in the package (confirmed by source grep) - standalone public API,
not dead/broken code.
[`validate_wkt_geometry`](https://jjmaynard.github.io/soilSIM/reference/validate_wkt_geometry.md)
(the package's live WKT validation entry point) does its own area check
inline rather than delegating here, since its area calculation
intentionally differs from this function's bounding-box approximation
for geographic coordinates (see
[`parse_wkt_geometry`](https://jjmaynard.github.io/soilSIM/reference/parse_wkt_geometry.md)'s
Usage note for the same distinction) and its `area_limits$max_action`
semantics (warn vs. error) predate this function.
