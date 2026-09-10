# Generated from utils.R (P1b split): WKT / geometry validation helpers.


#' Validate WKT Geometry (Generic)
#'
#' Comprehensive WKT geometry validation using terra package.
#' Supports multiple coordinate systems, validation contexts, and customizable thresholds.
#'
#' @param wkt_string WKT geometry string to validate
#' @param crs Coordinate reference system (default: "epsg:4326")
#' @param validation_context Context for validation ("geographic", "projected", "custom")
#' @param area_limits List with min/max area limits (units depend on CRS)
#' @param complexity_limits List with vertex and part limits
#' @param bounds_check List with coordinate bounds (xmin, xmax, ymin, ymax)
#' @param strict_mode Enable strict validation checks
#' @param return_geometry Whether to return the parsed geometry object
#'
#' @return Comprehensive geometry validation results
#'
#' @section Implementation note:
#' This is the package's single WKT validation entry point. It delegates to
#' the granular validators below (`validate_wkt_string()`,
#' `validate_geometry_validity()`, `validate_geometry_complexity()`,
#' `validate_geographic_context()`/`validate_projected_context()`) rather than
#' duplicating their logic inline. `validation_context`, `complexity_limits`,
#' and `strict_mode` are all read. The result is the nested
#' `geometry_stats$complexity_validation$complexity_stats` shape that
#' `adapter-ssurgo-acquire.R`'s downstream code expects.
#'
#' @export
validate_wkt_geometry <- function(wkt_string,
                                  crs = "epsg:4326",
                                  validation_context = "geographic",
                                  area_limits = list(min = 0.0001, max = 100, max_action = "warn"),
                                  complexity_limits = list(max_vertices = 1000, max_parts = 100),
                                  bounds_check = list(xmin = -180, xmax = 180, ymin = -90, ymax = 90),
                                  strict_mode = TRUE,
                                  return_geometry = FALSE) {

  log_message("DEBUG", paste("Starting WKT validation in", validation_context, "context"), category = "GeometryValidation")

  validation_result <- list(
    valid = TRUE,
    errors = character(0),
    warnings = character(0),
    geometry_stats = list()
  )

  # Basic string-format checks before attempting to parse - catches empty,
  # non-character, or unbalanced-parenthesis WKT with a clear message
  # instead of a raw terra/GDAL error. Reuses validate_wkt_string() rather
  # than duplicating this logic.
  string_check <- validate_wkt_string(wkt_string)
  validation_result$warnings <- c(validation_result$warnings, string_check$warnings)
  if (!string_check$valid) {
    validation_result$valid <- FALSE
    validation_result$errors <- c(validation_result$errors, string_check$errors)
    return(validation_result)
  }

  result <- tryCatch({
    log_message("DEBUG", paste("Parsing WKT geometry with CRS:", crs), category = "GeometryValidation")

    # Parse WKT with proper error handling
    if (!requireNamespace("terra", quietly = TRUE)) {
      stop("terra package required for geometry validation")
    }

    # Create geometry object
    geom <- terra::vect(wkt_string, crs = crs)

    # Check geometric validity - reuses validate_geometry_validity() so
    # `strict_mode` is honored consistently with the rest of this family.
    validity_check <- validate_geometry_validity(geom, strict_mode = strict_mode)
    validation_result$warnings <- c(validation_result$warnings, validity_check$warnings)
    if (!validity_check$valid) {
      validation_result$errors <- c(validation_result$errors, validity_check$errors)
      validation_result$valid <- FALSE
      return(validation_result)
    }

    # Extract coordinates safely
    coords_matrix <- tryCatch({
      terra::crds(geom)
    }, error = function(e) {
      log_message("DEBUG", paste("Coordinate extraction error:", e$message), category = "GeometryValidation")
      return(NULL)
    })

    if (is.null(coords_matrix) || !is.matrix(coords_matrix)) {
      validation_result$warnings <- c(validation_result$warnings, "Could not extract coordinates for detailed validation")
      validation_result$geometry_stats <- list(
        bbox = c(-180, -90, 180, 90),  # Default bounds
        area_deg2 = 1.0,  # Default area
        area_units = "degrees^2",
        geometry_type = "POLYGON",
        n_vertices = 4
      )
    } else {
      # Calculate bounding box
      bbox <- c(min(coords_matrix[,1]), min(coords_matrix[,2]),
                max(coords_matrix[,1]), max(coords_matrix[,2]))

      # Calculate area safely
      area_value <- tryCatch({
        as.numeric(terra::expanse(geom, unit = "m"))  # Get area in square meters
      }, error = function(e) {
        log_message("DEBUG", paste("Area calculation error:", e$message), category = "GeometryValidation")
        return(1.0)  # Default area
      })

      # Convert to degrees squared for geographic context
      if (is.finite(area_value) && area_value > 0) {
        # Rough conversion for small areas (not precise but reasonable for validation)
        area_deg2 <- area_value / (111000 * 111000)  # Approximate meters to degrees conversion
      } else {
        area_deg2 <- 1.0  # Default
      }

      # Store geometry statistics
      validation_result$geometry_stats <- list(
        bbox = bbox,
        area_deg2 = as.numeric(area_deg2),  # Ensure numeric
        area_value = as.numeric(area_value),  # Ensure numeric
        area_units = "degrees^2",
        geometry_type = as.character(terra::geomtype(geom)),
        n_vertices = nrow(coords_matrix)
      )
    }

    # Complexity check - passes `complexity_limits` through to validate_geometry_complexity().
    # Warning-only, so it can never flip an otherwise-valid geometry to invalid.
    complexity_check <- validate_geometry_complexity(geom, complexity_limits)
    validation_result$geometry_stats$complexity_validation <- complexity_check
    validation_result$warnings <- c(validation_result$warnings, complexity_check$warnings)

    # Context-specific checks - passes `validation_context` through to
    # validate_geographic_context()/validate_projected_context().
    context_check <- switch(validation_context,
      "geographic" = validate_geographic_context(geom, strict_mode = strict_mode),
      "projected" = validate_projected_context(geom, crs, strict_mode = strict_mode),
      NULL
    )
    if (!is.null(context_check)) {
      validation_result$warnings <- c(validation_result$warnings, context_check$warnings)
      if (!is.null(context_check$errors) && length(context_check$errors) > 0) {
        validation_result$errors <- c(validation_result$errors, context_check$errors)
        validation_result$valid <- FALSE
      }
    }

    # Bounds checking
    if (!is.null(validation_result$geometry_stats$bbox) && length(validation_result$geometry_stats$bbox) == 4) {
      bbox <- validation_result$geometry_stats$bbox
      if (bbox[1] < bounds_check$xmin || bbox[3] > bounds_check$xmax ||
          bbox[2] < bounds_check$ymin || bbox[4] > bounds_check$ymax) {
        validation_result$warnings <- c(validation_result$warnings, "Geometry extends beyond expected bounds")
      }
    }

    # Area validation
    if (!is.null(validation_result$geometry_stats$area_deg2) &&
        is.finite(validation_result$geometry_stats$area_deg2)) {
      area <- validation_result$geometry_stats$area_deg2
      if (area < area_limits$min) {
        validation_result$warnings <- c(validation_result$warnings, "Area below minimum threshold")
      } else if (area > area_limits$max) {
        if (area_limits$max_action == "warn") {
          validation_result$warnings <- c(validation_result$warnings, "Area above maximum threshold")
        } else {
          validation_result$errors <- c(validation_result$errors, "Area exceeds maximum allowed")
          validation_result$valid <- FALSE
        }
      }
    }

    if (return_geometry) {
      validation_result$geometry <- geom
    }

    validation_result

  }, error = function(e) {
    log_message("ERROR", paste("WKT parsing error:", e$message), category = "GeometryValidation")
    validation_result$valid <- FALSE
    validation_result$errors <- c(validation_result$errors, paste("Geometry validation error:", e$message))

    # Provide default geometry stats to prevent downstream errors
    validation_result$geometry_stats <- list(
      bbox = c(-180, -90, 180, 90),
      area_deg2 = 1.0,
      area_units = "degrees^2",
      geometry_type = "POLYGON",
      n_vertices = 4
    )
    validation_result
  })

  log_message("DEBUG", paste("WKT validation complete - Valid:", result$valid,
                             "Errors:", length(result$errors),
                             "Warnings:", length(result$warnings)),
              category = "GeometryValidation")

  return(result)
}

#' Validate WKT String Format
#'
#' Basic validation of WKT string format and content
#'
#' @param wkt_string WKT string to validate
#'
#' @return String validation results
#'
#' @export
validate_wkt_string <- function(wkt_string) {

  validation_results <- list(
    valid = TRUE,
    errors = character(0),
    warnings = character(0)
  )

  # Check for missing/empty string
  if (missing(wkt_string) || is.null(wkt_string) || wkt_string == "") {
    validation_results$valid <- FALSE
    validation_results$errors <- c(validation_results$errors, "WKT string cannot be empty or missing")
    return(validation_results)
  }

  # Check if it's a character string
  if (!is.character(wkt_string)) {
    validation_results$valid <- FALSE
    validation_results$errors <- c(validation_results$errors, "WKT must be a character string")
    return(validation_results)
  }

  # Check for multiple strings
  if (length(wkt_string) > 1) {
    validation_results$warnings <- c(validation_results$warnings,
                                     "Multiple WKT strings provided; only the first will be used")
    wkt_string <- wkt_string[1]
  }

  # Basic WKT format validation
  wkt_trimmed <- trimws(wkt_string)

  # Check for basic WKT geometry types
  valid_geometry_types <- c("POINT", "LINESTRING", "POLYGON", "MULTIPOINT",
                            "MULTILINESTRING", "MULTIPOLYGON", "GEOMETRYCOLLECTION")

  has_valid_type <- any(sapply(valid_geometry_types, function(type) {
    grepl(paste0("^", type, "\\s*\\("), wkt_trimmed, ignore.case = TRUE)
  }))

  if (!has_valid_type) {
    validation_results$warnings <- c(validation_results$warnings,
                                     "WKT string does not start with a recognized geometry type")
  }

  # Check for balanced parentheses
  open_parens <- nchar(wkt_trimmed) - nchar(gsub("\\(", "", wkt_trimmed))
  close_parens <- nchar(wkt_trimmed) - nchar(gsub("\\)", "", wkt_trimmed))

  if (open_parens != close_parens) {
    validation_results$valid <- FALSE
    validation_results$errors <- c(validation_results$errors,
                                   paste("Unbalanced parentheses in WKT string (", open_parens, "open,", close_parens, "close)"))
  }

  # Check string length (very long WKT strings may cause performance issues)
  if (nchar(wkt_string) > 100000) {
    validation_results$warnings <- c(validation_results$warnings,
                                     "Very long WKT string may cause performance issues")
  }

  return(validation_results)
}

#' Parse WKT Geometry
#'
#' Parse WKT string into terra geometry object with error handling
#'
#' @param wkt_string Valid WKT string
#' @param crs Coordinate reference system
#'
#' @return Parsing results with geometry object and statistics
#'
#' @section Usage note:
#' Fully implemented and exported, but not currently called from anywhere
#' else in the package (confirmed by source grep) - standalone public API,
#' not dead/broken code. \code{\link{validate_wkt_geometry}} (the package's
#' live WKT validation entry point) parses independently rather than
#' delegating here, since its area-statistics calculation intentionally
#' differs from this function's bounding-box approximation for geographic
#' coordinates (it always computes true polygon area via
#' \code{terra::expanse()}, then converts to degrees^2 - see that
#' function's own comments).
#'
#' @export
parse_wkt_geometry <- function(wkt_string, crs) {

  parsing_results <- list(
    valid = TRUE,
    errors = character(0),
    warnings = character(0),
    geometry = NULL,
    stats = list()
  )

  log_message("DEBUG", paste("Parsing WKT geometry with CRS:", crs), category = "GeometryValidation")

  tryCatch({
    # Parse WKT with terra
    geom <- terra::vect(wkt_string, crs = crs)

    parsing_results$geometry <- geom

    # Calculate basic geometry statistics
    bbox <- terra::ext(geom)

    parsing_results$stats <- list(
      geometry_type = terra::geomtype(geom),
      n_geometries = terra::nrow(geom),
      n_vertices = sum(terra::nrow(terra::geom(geom))),
      bbox = as.numeric(bbox),
      crs = terra::crs(geom),
      is_valid = terra::is.valid(geom)
    )

    # Calculate area (handling different CRS units)
    if (terra::is.lonlat(geom)) {
      # Geographic coordinates - area in square degrees (approximate)
      area_value <- (bbox[2] - bbox[1]) * (bbox[4] - bbox[3])
      area_units <- "deg^2"
    } else {
      # Projected coordinates - use terra's area calculation
      area_value <- terra::expanse(geom, unit = "m")
      area_units <- "m^2"
    }

    parsing_results$stats$area_value <- area_value
    parsing_results$stats$area_units <- area_units

    log_message("DEBUG", paste("Geometry parsed successfully - Type:",
                               parsing_results$stats$geometry_type,
                               ", Vertices:", parsing_results$stats$n_vertices,
                               ", Area:", round(area_value, 4), area_units),
                category = "GeometryValidation")

  }, error = function(e) {
    parsing_results$valid <- FALSE
    parsing_results$errors <- c(parsing_results$errors,
                                paste("WKT parsing failed:", e$message))

    log_message("ERROR", paste("WKT parsing error:", e$message), category = "GeometryValidation")

    # Try to provide more specific error information
    if (grepl("GDAL", e$message, ignore.case = TRUE)) {
      parsing_results$errors <- c(parsing_results$errors,
                                  "GDAL error - check WKT syntax and coordinate values")
    }
    if (grepl("CRS", e$message, ignore.case = TRUE)) {
      parsing_results$warnings <- c(parsing_results$warnings,
                                    paste("CRS issue with", crs, "- consider using a different CRS"))
    }
  })

  return(parsing_results)
}

#' Validate Geometry Validity
#'
#' Check geometric validity using terra's validation
#'
#' @param geom Terra geometry object
#' @param strict_mode Enable strict validation
#'
#' @return Validity validation results
#'
#' @export
validate_geometry_validity <- function(geom, strict_mode = TRUE) {

  validity_results <- list(
    valid = TRUE,
    errors = character(0),
    warnings = character(0)
  )

  # Basic validity check
  if (!terra::is.valid(geom)) {
    if (strict_mode) {
      validity_results$valid <- FALSE
      validity_results$errors <- c(validity_results$errors, "Geometry is not topologically valid")
    } else {
      validity_results$warnings <- c(validity_results$warnings, "Geometry has validity issues but proceeding")
    }

    # Try to get more specific validity information
    tryCatch({
      # Attempt to fix and compare
      fixed_geom <- terra::makeValid(geom)
      if (terra::is.valid(fixed_geom)) {
        validity_results$warnings <- c(validity_results$warnings,
                                       "Geometry can be automatically fixed with terra::makeValid()")
      }
    }, error = function(e) {
      validity_results$warnings <- c(validity_results$warnings,
                                     "Cannot automatically fix geometry validity issues")
    })
  }

  return(validity_results)
}

#' Validate Geometry Area
#'
#' Validate geometry area against specified limits
#'
#' @param geom Terra geometry object
#' @param area_limits List with min/max area limits
#' @param context Validation context for appropriate messaging
#'
#' @return Area validation results
#'
#' @section Usage note:
#' Fully implemented and exported, but not currently called from anywhere
#' else in the package (confirmed by source grep) - standalone public API,
#' not dead/broken code. \code{\link{validate_wkt_geometry}} (the package's
#' live WKT validation entry point) does its own area check inline rather
#' than delegating here, since its area calculation intentionally differs
#' from this function's bounding-box approximation for geographic
#' coordinates (see \code{\link{parse_wkt_geometry}}'s Usage note for the
#' same distinction) and its `area_limits$max_action` semantics (warn vs.
#' error) predate this function.
#'
#' @export
validate_geometry_area <- function(geom, area_limits, context = "general") {

  area_results <- list(
    valid = TRUE,
    errors = character(0),
    warnings = character(0),
    area_stats = list()
  )

  # Calculate area
  if (terra::is.lonlat(geom)) {
    bbox <- terra::ext(geom)
    area_value <- (bbox[2] - bbox[1]) * (bbox[4] - bbox[3])
    area_units <- "deg^2"
  } else {
    area_value <- terra::expanse(geom, unit = "m")
    area_units <- "m^2"
  }

  area_results$area_stats <- list(
    area_value = area_value,
    area_units = area_units,
    area_limits = area_limits
  )

  # Validate against limits
  if (!is.null(area_limits$min) && area_value < area_limits$min) {
    if (!is.null(area_limits$min_action) && area_limits$min_action == "error") {
      area_results$valid <- FALSE
      area_results$errors <- c(area_results$errors,
                               paste("Geometry area too small:", round(area_value, 6), area_units,
                                     "< minimum", area_limits$min, area_units))
    } else {
      area_results$warnings <- c(area_results$warnings,
                                 paste("Geometry area is very small:", round(area_value, 6), area_units))
    }
  }

  if (!is.null(area_limits$max) && area_value > area_limits$max) {
    if (!is.null(area_limits$max_action) && area_limits$max_action == "error") {
      area_results$valid <- FALSE
      area_results$errors <- c(area_results$errors,
                               paste("Geometry area too large:", round(area_value, 2), area_units,
                                     "> maximum", area_limits$max, area_units))
    } else {
      area_results$warnings <- c(area_results$warnings,
                                 paste("Geometry area is very large:", round(area_value, 2), area_units,
                                       "- processing may be slow"))
    }
  }

  return(area_results)
}

#' Validate Geometry Complexity
#'
#' Check geometry complexity (vertices, parts, etc.)
#'
#' @param geom Terra geometry object
#' @param complexity_limits List with complexity thresholds
#'
#' @return Complexity validation results
#'
#' @export
validate_geometry_complexity <- function(geom, complexity_limits) {

  complexity_results <- list(
    warnings = character(0),
    complexity_stats = list()
  )

  # Calculate complexity metrics
  n_vertices <- sum(terra::nrow(terra::geom(geom)))
  n_parts <- terra::nrow(geom)
  geometry_type <- terra::geomtype(geom)

  complexity_results$complexity_stats <- list(
    n_vertices = n_vertices,
    n_parts = n_parts,
    geometry_type = geometry_type,
    complexity_score = calculate_complexity_score(n_vertices, n_parts, geometry_type)
  )

  # Check against limits
  if (!is.null(complexity_limits$max_vertices) && n_vertices > complexity_limits$max_vertices) {
    complexity_results$warnings <- c(complexity_results$warnings,
                                     paste("Geometry has many vertices (", n_vertices, ") - consider simplifying for better performance"))
  }

  if (!is.null(complexity_limits$max_parts) && n_parts > complexity_limits$max_parts) {
    complexity_results$warnings <- c(complexity_results$warnings,
                                     paste("Geometry has many parts (", n_parts, ") - processing may be slow"))
  }

  return(complexity_results)
}

#' Validate Coordinate Bounds
#'
#' Check if coordinates fall within expected bounds
#'
#' @param geom Terra geometry object
#' @param bounds_check List with coordinate bounds
#' @param context Validation context
#'
#' @return Bounds validation results
#'
#' @section Usage note:
#' Fully implemented and exported, but not currently called from anywhere
#' else in the package (confirmed by source grep) - standalone public API,
#' not dead/broken code. \code{\link{validate_wkt_geometry}} (the package's
#' live WKT validation entry point) does its own bounds check inline rather
#' than delegating here: this function treats an out-of-bounds geometry as
#' an error (`valid = FALSE`), while the live path treats it as a
#' warning-only condition - adopting this function's stricter behavior
#' there would be a real behavior change for existing callers, not a pure
#' cleanup.
#'
#' @export
validate_coordinate_bounds <- function(geom, bounds_check, context = "general") {

  bounds_results <- list(
    valid = TRUE,
    errors = character(0),
    warnings = character(0)
  )

  bbox <- terra::ext(geom)

  # Check X bounds
  if (bbox[1] < bounds_check$xmin || bbox[2] > bounds_check$xmax) {
    bounds_results$valid <- FALSE
    bounds_results$errors <- c(bounds_results$errors,
                               paste("X coordinates outside valid bounds:",
                                     round(bbox[1], 4), "to", round(bbox[2], 4),
                                     "not within", bounds_check$xmin, "to", bounds_check$xmax))
  }

  # Check Y bounds
  if (bbox[3] < bounds_check$ymin || bbox[4] > bounds_check$ymax) {
    bounds_results$valid <- FALSE
    bounds_results$errors <- c(bounds_results$errors,
                               paste("Y coordinates outside valid bounds:",
                                     round(bbox[3], 4), "to", round(bbox[4], 4),
                                     "not within", bounds_check$ymin, "to", bounds_check$ymax))
  }

  return(bounds_results)
}

#' Get Validation Defaults by Context
#'
#' Get appropriate default validation parameters for different contexts
#'
#' @param context Validation context
#' @param crs Coordinate reference system
#'
#' @return Default validation parameters
#'
#' @section Usage note:
#' Fully implemented and exported, but not currently called from anywhere
#' else in the package (confirmed by source grep) - standalone public API,
#' not dead/broken code. \code{\link{validate_wkt_geometry}} (the package's
#' live WKT validation entry point) hardcodes its own context-appropriate
#' defaults directly in its function signature rather than calling this at
#' runtime, so the two default sets should be kept in sync by hand if either
#' changes.
#'
#' @export
get_validation_defaults <- function(context, crs) {

  defaults <- switch(context,
                     "geographic" = list(
                       area_limits = list(min = 0.0001, max = 100, max_action = "warn"),
                       complexity_limits = list(max_vertices = 1000, max_parts = 100),
                       bounds_check = list(xmin = -180, xmax = 180, ymin = -90, ymax = 90)
                     ),

                     "projected" = list(
                       area_limits = list(min = 1, max = 1e12, max_action = "warn"), # 1 m^2 to 1M km^2
                       complexity_limits = list(max_vertices = 10000, max_parts = 1000),
                       bounds_check = NULL  # Depends on specific projection
                     ),

                     "custom" = list(
                       area_limits = NULL,
                       complexity_limits = list(max_vertices = 5000, max_parts = 500),
                       bounds_check = NULL
                     ),

                     # Default
                     list(
                       area_limits = NULL,
                       complexity_limits = list(max_vertices = 1000, max_parts = 100),
                       bounds_check = NULL
                     )
  )

  return(defaults)
}

#' Calculate Complexity Score
#'
#' Calculate a complexity score for geometry
#'
#' @param n_vertices Number of vertices
#' @param n_parts Number of parts/geometries
#' @param geometry_type Type of geometry
#'
#' @return Complexity score
#'
#' @export
calculate_complexity_score <- function(n_vertices, n_parts, geometry_type) {

  # Base score from vertices and parts
  vertex_score <- log10(max(n_vertices, 1))
  part_score <- log10(max(n_parts, 1))

  # Type complexity multiplier
  type_multiplier <- switch(geometry_type,
                            "POINT" = 1,
                            "LINESTRING" = 1.2,
                            "POLYGON" = 1.5,
                            "MULTIPOINT" = 1.3,
                            "MULTILINESTRING" = 1.6,
                            "MULTIPOLYGON" = 2.0,
                            "GEOMETRYCOLLECTION" = 2.5,
                            1.0
  )

  complexity_score <- (vertex_score + part_score) * type_multiplier

  return(round(complexity_score, 2))
}

#' Validate Geographic Context
#'
#' Additional validations specific to geographic (lat/lon) coordinates
#'
#' @param geom Terra geometry object
#' @param strict_mode Enable strict validation
#'
#' @return Geographic validation results
#'
#' @export
validate_geographic_context <- function(geom, strict_mode = TRUE) {

  geo_results <- list(
    errors = character(0),
    warnings = character(0)
  )

  if (!terra::is.lonlat(geom)) {
    geo_results$warnings <- c(geo_results$warnings,
                              "Geometry is not in geographic coordinates but geographic context was specified")
    return(geo_results)
  }

  bbox <- terra::ext(geom)

  # Check for impossible geographic coordinates
  if (strict_mode) {
    if (abs(bbox[1]) > 180 || abs(bbox[2]) > 180) {
      geo_results$errors <- c(geo_results$errors,
                              "Longitude values outside valid range (-180 to 180)")
    }

    if (abs(bbox[3]) > 90 || abs(bbox[4]) > 90) {
      geo_results$errors <- c(geo_results$errors,
                              "Latitude values outside valid range (-90 to 90)")
    }
  }

  # Check for suspicious coordinate patterns
  if (bbox[1] == bbox[2] && bbox[3] == bbox[4]) {
    geo_results$warnings <- c(geo_results$warnings,
                              "Geometry appears to be a single point")
  }

  return(geo_results)
}

#' Validate Projected Context
#'
#' Additional validations for projected coordinate systems
#'
#' @param geom Terra geometry object
#' @param crs Coordinate reference system
#' @param strict_mode Enable strict validation
#'
#' @return Projected validation results
#'
#' @export
validate_projected_context <- function(geom, crs, strict_mode = TRUE) {

  proj_results <- list(
    warnings = character(0)
  )

  if (terra::is.lonlat(geom)) {
    proj_results$warnings <- c(proj_results$warnings,
                               "Geometry is in geographic coordinates but projected context was specified")
  }

  # Check for very large coordinate values that might indicate incorrect CRS
  bbox <- terra::ext(geom)
  max_coord <- max(abs(c(bbox[1], bbox[2], bbox[3], bbox[4])))

  if (max_coord > 1e8) {
    proj_results$warnings <- c(proj_results$warnings,
                               "Very large coordinate values detected - verify CRS is correct")
  }

  return(proj_results)
}
