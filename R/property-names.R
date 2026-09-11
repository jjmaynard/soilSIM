# Generated from utils.R (P1b split): soil-property name validation, synonyms & lookup


#' Validate Soil Properties (Simple Generic Version)
#'
#' Simple, extensible property validation that works with any property source.
#' Focuses on core validation logic without excessive complexity.
#'
#' @param properties Vector of property names to validate
#' @param property_lookup Property lookup source (function, vector, data.frame, or source name)
#' @param strict_mode Enable strict validation (default: TRUE)
#' @param max_invalid_pct Maximum percentage of invalid properties allowed (default: 50)
#' @param performance_threshold Warn if more than this many properties (default: 10)
#'
#' @return Simple validation results list
#'
#' @keywords internal
validate_properties <- function(properties,
                                property_lookup = "ssurgo",
                                strict_mode = TRUE,
                                max_invalid_pct = 50,
                                performance_threshold = 10) {

  log_message("DEBUG", paste("Validating", length(properties), "properties"), category = "PropertyValidation")

  # Initialize results
  results <- list(
    valid = TRUE,
    errors = character(0),
    warnings = character(0),
    property_stats = list()
  )

  # 1. Basic input validation
  if (!is.character(properties) || length(properties) == 0) {
    results$valid <- FALSE
    results$errors <- c(results$errors, "Properties must be a non-empty character vector")
    return(results)
  }

  # 2. Get available properties
  available_props <- available_properties(property_lookup)

  if (length(available_props) == 0) {
    results$warnings <- c(results$warnings, "Could not load property lookup - skipping validation")
    results$property_stats$validation_method <- "skipped"
    return(results)
  }

  # 3. Match properties
  valid_props <- intersect(properties, available_props)
  invalid_props <- setdiff(properties, available_props)

  # 4. Calculate statistics
  results$property_stats <- list(
    total_requested = length(properties),
    valid_properties = length(valid_props),
    invalid_properties = length(invalid_props),
    validation_method = "direct_match",
    invalid_property_names = if(length(invalid_props) > 0) invalid_props else NULL
  )

  # 5. Apply validation rules
  if (length(invalid_props) > 0) {
    results$warnings <- c(results$warnings,
                          paste("Unknown properties will be ignored:", paste(invalid_props, collapse = ", ")))
  }

  # Strict mode: fail if too many invalid properties
  invalid_pct <- (length(invalid_props) / length(properties)) * 100
  if (strict_mode && invalid_pct > max_invalid_pct) {
    results$valid <- FALSE
    results$errors <- c(results$errors,
                        paste("Too many invalid properties:", length(invalid_props), "of", length(properties),
                              paste0("(", round(invalid_pct, 1), "% > ", max_invalid_pct, "%)")))
  }

  # Performance warning
  if (length(valid_props) > performance_threshold) {
    results$warnings <- c(results$warnings,
                          paste("Large number of properties (", length(valid_props), ") may impact performance"))
  }

  # No valid properties
  if (length(valid_props) == 0) {
    results$valid <- FALSE
    results$errors <- c(results$errors, "No valid properties found")
  }

  log_message("DEBUG", paste("Property validation complete - Valid:", length(valid_props),
                             "Invalid:", length(invalid_props)), category = "PropertyValidation")

  return(results)
}

#' Get Available Properties from Various Sources
#'
#' Simple property lookup that handles multiple source types
#'
#' @param property_lookup Property source (string, function, vector, or data.frame)
#'
#' @return Vector of available property names
#'
#' @export
available_properties <- function(property_lookup) {

  tryCatch({

    if (is.character(property_lookup) && length(property_lookup) == 1) {
      # Named source
      return(predefined_properties(property_lookup))

    } else if (is.function(property_lookup)) {
      # Function that returns properties
      return(as.character(property_lookup()))

    } else if (is.character(property_lookup)) {
      # Vector of property names
      return(property_lookup)

    } else if (is.data.frame(property_lookup)) {
      # Data frame with properties
      return(extract_properties_from_dataframe(property_lookup))

    } else {
      log_message("WARN", "Unsupported property lookup type", category = "PropertyValidation")
      return(character(0))
    }

  }, error = function(e) {
    log_message("WARN", paste("Error loading properties:", e$message), category = "PropertyValidation")
    return(character(0))
  })
}

#' Get Predefined Property Sets
#'
#' Returns property lists for common soil data sources
#'
#' @param source_name Name of the property source
#'
#' @return Vector of property names
#'
#' @export
predefined_properties <- function(source_name) {

  switch(tolower(source_name),

         "ssurgo" = {
           # build_ssurgo_property_lookup() lives in adapter-ssurgo-acquire.R. The
           # tryCatch() below degrades to a WARN plus an empty vector, since this switch
           # has no other way to signal an unknown or unavailable source.
           tryCatch({
             ssurgo_data <- build_ssurgo_property_lookup()
             if ("Property" %in% names(ssurgo_data)) {
               return(ssurgo_data$Property)
             } else {
               return(character(0))
             }
           }, error = function(e) {
             log_message("WARN", "Could not load SSURGO properties", category = "PropertyValidation")
             return(character(0))
           })
         },

         "nrcs" = c(
           "clay_pct", "sand_pct", "silt_pct", "bulk_density", "ph",
           "organic_matter", "cec", "available_water_capacity"
         ),

         "laboratory" = c(
           "clay", "sand", "silt", "bulk_density", "ph_water", "ph_cacl2",
           "organic_carbon", "total_nitrogen", "cec", "base_saturation"
         ),

         "basic" = c(
           "clay", "sand", "silt", "bulk_density", "ph", "organic_matter"
         ),

         # Default: return empty vector
         {
           log_message("WARN", paste("Unknown property source:", source_name), category = "PropertyValidation")
           character(0)
         }
  )
}

#' Extract Properties from Data Frame
#'
#' Extracts property names from a data frame using common column names
#'
#' @param df Data frame containing property information
#'
#' @return Vector of property names
#'
#' @keywords internal
extract_properties_from_dataframe <- function(df) {

  # Try common column names for properties
  property_columns <- c("property", "Property", "name", "variable", "column")

  for (col in property_columns) {
    if (col %in% names(df)) {
      return(as.character(df[[col]]))
    }
  }

  # Fallback: use first column
  log_message("DEBUG", "Using first column as property names", category = "PropertyValidation")
  return(as.character(df[[1]]))
}

#' Validate Properties with Synonyms
#'
#'
#' @param properties Vector of property names to validate
#' @param property_lookup Property lookup source
#' @param strict_mode Enable strict validation
#' @param max_invalid_pct Maximum percentage of invalid properties allowed
#' @param performance_threshold Performance warning threshold
#'
#' @return validation results with synonym matching
#'
#' @keywords internal
validate_properties_with_synonyms <- function(properties,
                                              property_lookup = "ssurgo",
                                              strict_mode = TRUE,
                                              max_invalid_pct = 50,
                                              performance_threshold = 10) {

  log_message("DEBUG", "Validating properties with synonym matching", category = "PropertyValidation")
  log_message("DEBUG", paste("Validating", length(properties), "properties"), category = "PropertyValidation")

  validation_result <- list(
    valid = TRUE,
    errors = character(0),
    warnings = character(0),
    property_stats = list()
  )

  # Known SSURGO properties and their synonyms, derived from the shared
  # canonical table (see .ssurgo_property_synonyms()) so this list can't
  # drift out of sync with get_default_property_mapping()/default_property_synonyms().
  ssurgo_canonical <- .ssurgo_property_synonyms()
  ssurgo_properties <- stats::setNames(
    lapply(names(ssurgo_canonical), function(canonical_name) {
      base_name <- sub("_r$", "", canonical_name)
      unique(c(base_name, canonical_name, ssurgo_canonical[[canonical_name]]))
    }),
    sub("_r$", "", names(ssurgo_canonical))
  )

  # Validate against known properties
  valid_properties <- character(0)
  invalid_properties <- character(0)
  synonym_matches <- 0

  for (prop in properties) {
    # Check direct matches first
    direct_match <- FALSE
    for (known_prop in names(ssurgo_properties)) {
      if (prop %in% ssurgo_properties[[known_prop]]) {
        valid_properties <- c(valid_properties, prop)
        direct_match <- TRUE
        if (prop != known_prop) {
          synonym_matches <- synonym_matches + 1
        }
        break
      }
    }

    if (!direct_match) {
      # Check for partial matches (e.g., property names with suffixes)
      base_prop <- gsub("_(r|l|h)$", "", prop)  # Remove _r, _l, _h suffixes

      partial_match <- FALSE
      for (known_prop in names(ssurgo_properties)) {
        if (base_prop %in% ssurgo_properties[[known_prop]] ||
            known_prop == base_prop) {
          valid_properties <- c(valid_properties, prop)
          partial_match <- TRUE
          synonym_matches <- synonym_matches + 1
          break
        }
      }

      if (!partial_match) {
        invalid_properties <- c(invalid_properties, prop)
      }
    }
  }

  # Calculate validation statistics
  total_properties <- length(properties)
  valid_count <- length(valid_properties)
  invalid_count <- length(invalid_properties)
  success_rate <- valid_count / total_properties

  log_message("DEBUG", paste("Property validation complete - Valid:", valid_count, "Invalid:", invalid_count),
              category = "PropertyValidation")
  log_message("DEBUG", paste("Matched", synonym_matches, "properties using synonyms"),
              category = "PropertyValidation")

  # Store validation statistics
  validation_result$property_stats <- list(
    total_requested = total_properties,
    valid_properties = valid_count,
    invalid_properties = invalid_count,
    success_rate = success_rate,
    synonym_matches = synonym_matches,
    validation_method = "enhanced_with_synonyms",
    valid_property_names = valid_properties,
    invalid_property_names = invalid_properties
  )

  # Determine validation result
  invalid_percentage <- (invalid_count / total_properties) * 100

  if (invalid_count > 0) {
    if (invalid_percentage > max_invalid_pct) {
      validation_result$valid <- FALSE
      validation_result$errors <- c(validation_result$errors,
                                    paste("Too many invalid properties:", invalid_count, "of", total_properties))
    } else {
      validation_result$warnings <- c(validation_result$warnings,
                                      paste("Some invalid properties:", paste(invalid_properties, collapse = ", ")))
    }
  }

  # Performance warning
  if (total_properties > performance_threshold) {
    validation_result$warnings <- c(validation_result$warnings,
                                    paste("Large number of properties may impact performance:", total_properties))
  }

  return(validation_result)
}

#' Get Default Synonyms for Property Sources
#'
#' Returns common synonym mappings for different property sources
#'
#' @param property_lookup Property lookup source
#'
#' @return Named list of synonyms
#'
#' @section Usage note:
#' Fully implemented and exported, but not currently called from anywhere
#' else in the package (confirmed by source grep) - standalone public API
#' for callers who need it, not dead/broken code. Its `"ssurgo"` branch
#' shares the same canonical synonym data as \code{get_default_property_mapping()}
#' (internal) and \code{\link{validate_properties_with_synonyms}} (see
#' `.ssurgo_property_synonyms()`).
#'
#' @export
default_property_synonyms <- function(property_lookup) {

  if (is.character(property_lookup) && length(property_lookup) == 1) {
    source_name <- tolower(property_lookup)
  } else {
    source_name <- "general"
  }

  switch(source_name,

         "ssurgo" = .ssurgo_property_synonyms()[c("claytotal_r", "sandtotal_r", "silttotal_r",
                                                   "dbovendry_r", "ph1to1h2o_r", "om_r", "cec7_r")],

         "nrcs" = list(
           "clay_pct" = c("clay", "claytotal_r"),
           "sand_pct" = c("sand", "sandtotal_r"),
           "silt_pct" = c("silt", "silttotal_r"),
           "bulk_density" = c("bd", "dbovendry_r"),
           "ph" = c("pH", "ph1to1h2o_r", "ph_water"),
           "organic_matter" = c("om", "om_r"),
           "cec" = c("cec7_r", "cation_exchange_capacity")
         ),

         "laboratory" = list(
           "clay" = c("clay_pct", "particle_size_clay"),
           "sand" = c("sand_pct", "particle_size_sand"),
           "silt" = c("silt_pct", "particle_size_silt"),
           "bulk_density" = c("bd", "bulk_density_field"),
           "ph_water" = c("ph", "pH", "ph_h2o"),
           "organic_carbon" = c("oc", "soc", "carbon_organic")
         ),

         # Default: empty list
         list()
  )
}

#' Simple Property Lookup Creator
#'
#' Helper function to create custom property lookups easily
#'
#' @param properties Vector of property names
#' @param synonyms Optional named list of synonyms
#' @param metadata Optional metadata list
#'
#' @return Property lookup object
#'
#' @section Usage note:
#' Fully implemented and exported, but not currently called from anywhere
#' else in the package (confirmed by source grep) - standalone public API
#' for callers who need it, not dead/broken code.
#'
#' @keywords internal
create_property_lookup <- function(properties, synonyms = NULL, metadata = NULL) {

  lookup <- list(
    properties = as.character(properties),
    synonyms = synonyms,
    metadata = metadata
  )

  # Add a function interface
  lookup$get_properties <- function() {
    return(lookup$properties)
  }

  lookup$validate <- function(props, strict_mode = TRUE) {
    if (is.null(lookup$synonyms)) {
      return(validate_properties(props, lookup$properties, strict_mode))
    } else {
      # validate_properties_with_synonyms() has no `synonyms` parameter (it uses its own
      # internal SSURGO synonym table). Pass strict_mode by name to match the actual signature.
      return(validate_properties_with_synonyms(props, lookup$properties, strict_mode = strict_mode))
    }
  }

  class(lookup) <- "property_lookup"
  return(lookup)
}

