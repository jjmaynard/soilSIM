# Validate Coordinate Bounds

Check if coordinates fall within expected bounds

## Usage

``` r
validate_coordinate_bounds(geom, bounds_check, context = "general")
```

## Arguments

- geom:

  Terra geometry object

- bounds_check:

  List with coordinate bounds

- context:

  Validation context

## Value

Bounds validation results

## Usage note

Fully implemented and exported, but not currently called from anywhere
else in the package (confirmed by source grep) - standalone public API,
not dead/broken code.
[`validate_wkt_geometry`](https://jjmaynard.github.io/soilSIM/reference/validate_wkt_geometry.md)
(the package's live WKT validation entry point) does its own bounds
check inline rather than delegating here: this function treats an
out-of-bounds geometry as an error (`valid = FALSE`), while the live
path treats it as a warning-only condition - adopting this function's
stricter behavior there would be a real behavior change for existing
callers, not a pure cleanup.
