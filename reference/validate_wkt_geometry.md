# Validate WKT Geometry (Generic)

Comprehensive WKT geometry validation using terra package. Supports
multiple coordinate systems, validation contexts, and customizable
thresholds.

## Usage

``` r
validate_wkt_geometry(
  wkt_string,
  crs = "epsg:4326",
  validation_context = "geographic",
  area_limits = list(min = 1e-04, max = 100, max_action = "warn"),
  complexity_limits = list(max_vertices = 1000, max_parts = 100),
  bounds_check = list(xmin = -180, xmax = 180, ymin = -90, ymax = 90),
  strict_mode = TRUE,
  return_geometry = FALSE
)
```

## Arguments

- wkt_string:

  WKT geometry string to validate

- crs:

  Coordinate reference system (default: "epsg:4326")

- validation_context:

  Context for validation ("geographic", "projected", "custom")

- area_limits:

  List with min/max area limits (units depend on CRS)

- complexity_limits:

  List with vertex and part limits

- bounds_check:

  List with coordinate bounds (xmin, xmax, ymin, ymax)

- strict_mode:

  Enable strict validation checks

- return_geometry:

  Whether to return the parsed geometry object

## Value

Comprehensive geometry validation results

## Implementation note

This is the package's single WKT validation entry point. It delegates to
the granular validators below
([`validate_wkt_string()`](https://jjmaynard.github.io/soilSIM/reference/validate_wkt_string.md),
[`validate_geometry_validity()`](https://jjmaynard.github.io/soilSIM/reference/validate_geometry_validity.md),
[`validate_geometry_complexity()`](https://jjmaynard.github.io/soilSIM/reference/validate_geometry_complexity.md),
[`validate_geographic_context()`](https://jjmaynard.github.io/soilSIM/reference/validate_geographic_context.md)/[`validate_projected_context()`](https://jjmaynard.github.io/soilSIM/reference/validate_projected_context.md))
rather than duplicating their logic inline. `validation_context`,
`complexity_limits`, and `strict_mode` are all read. The result is the
nested `geometry_stats$complexity_validation$complexity_stats` shape
that `adapter-ssurgo-acquire.R`'s downstream code expects.
