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
