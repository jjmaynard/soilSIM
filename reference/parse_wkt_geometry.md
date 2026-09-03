# Parse WKT Geometry

Parse WKT string into terra geometry object with error handling

## Usage

``` r
parse_wkt_geometry(wkt_string, crs)
```

## Arguments

- wkt_string:

  Valid WKT string

- crs:

  Coordinate reference system

## Value

Parsing results with geometry object and statistics

## Usage note

Fully implemented and exported, but not currently called from anywhere
else in the package (confirmed by source grep) - standalone public API,
not dead/broken code.
[`validate_wkt_geometry`](https://jjmaynard.github.io/soilSIM/reference/validate_wkt_geometry.md)
(the package's live WKT validation entry point) parses independently
rather than delegating here, since its area-statistics calculation
intentionally differs from this function's bounding-box approximation
for geographic coordinates (it always computes true polygon area via
[`terra::expanse()`](https://rspatial.github.io/terra/reference/expanse.html),
then converts to degrees^2 - see that function's own comments).
