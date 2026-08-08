#' @title Soil Property Infilling System
#' @description Comprehensive soil property infilling with automatic exclusion of unsuitable horizons
#' @name soil_infilling
NULL

# ==============================================================================
# 1. MAIN INFILLING FUNCTIONS
# ==============================================================================

#' Comprehensive Soil Property Processing
#'
#' Complete soil property processing workflow that handles multiple properties
#' with automatic exclusion of unsuitable horizons and intelligent infilling strategies.
#'
#' @param df Input soil data frame
#' @param properties Vector of properties to process (NULL = auto-detect)
#' @param max_depth Maximum depth for processing (default: 250 cm)
#' @param remove_unsuitable Whether to remove unsuitable horizons from output
#' @param remove_incomplete Whether to remove incomplete rows
#' @param required_properties Vector of properties that must be complete
#' @param verbose Whether to print detailed progress messages
#'
#' @return Data frame with processed and infilled properties
#'
#' @export
process_soil_properties_comprehensive <- function(df,
                                                  properties = NULL,
                                                  max_depth = 250,
                                                  remove_unsuitable = FALSE,
                                                  remove_incomplete = FALSE,
                                                  required_properties = NULL,
                                                  verbose = getOption("ssurgo.verbose", FALSE)) {
  .old_log_cfg <- set_verbose_logging(verbose)
  on.exit(options(soil_workflow_log_config = .old_log_cfg), add = TRUE)

  start_time <- Sys.time()

  if (verbose) {
    log_message("INFO", "=== COMPREHENSIVE SOIL PROPERTY PROCESSING ===", category = "Infilling")
  }

  # Step 1: Auto-detect properties if not specified
  if (is.null(properties)) {
    properties <- auto_detect_soil_properties(df)
    if (verbose) {
      log_message("INFO", paste("Auto-detected properties:", paste(properties, collapse = ", ")), category = "Infilling")
    }
  }

  # Validate properties using Module 0
  property_validation <- validate_properties_with_synonyms(
    properties,
    property_lookup = "ssurgo",
    strict_mode = FALSE
  )

  if (!property_validation$valid) {
    log_message("WARN", paste("Property validation warnings:",
                              paste(property_validation$warnings, collapse = "; ")), category = "Infilling")
  }

  # Step 2: Initialize tracking and ensure essential columns
  df <- ensure_infilling_columns(df, verbose = verbose)

  # Step 3: Process properties by category for optimal results
  processing_log <- list()

  # Phase 1: Foundation properties (texture, bulk density, rock fragments)
  foundation_props <- intersect(c("sandtotal", "claytotal", "silttotal", "dbovendry", "rfv"), properties)

  if (length(foundation_props) > 0) {
    if (verbose) log_message("INFO", "Phase 1: Processing foundation properties", category = "Infilling")

    for (prop in foundation_props) {
      if (verbose) log_message("DEBUG", paste("Processing", prop), category = "Infilling")

      result <- process_single_property(df, prop, max_depth, verbose)
      df <- result$data
      processing_log[[prop]] <- result$log
    }
  }

  # Phase 2: Water retention estimation (if texture + BD available)
  water_ret_props <- intersect(c("wthirdbar", "wfifteenbar"), properties)

  if (length(water_ret_props) > 0) {
    texture_bd_available <- all(paste0(c("sandtotal", "claytotal", "silttotal", "dbovendry"), "_r") %in% names(df))

    if (texture_bd_available) {
      if (verbose) log_message("INFO", "Phase 2: Estimating water retention using Saxton-Rawls", category = "Infilling")

      df <- infill_water_retention_saxton_rawls_integrated(df, max_depth = max_depth, verbose = verbose)

      # Track water retention results
      for (prop in water_ret_props) {
        r_col <- paste0(prop, "_r")
        if (r_col %in% names(df)) {
          processing_log[[prop]] <- list(
            method = "saxton_rawls",
            estimated_values = sum(!is.na(df[[r_col]]))
          )
        }
      }
    }
  }

  # Phase 3: Chemical properties
  chem_props <- setdiff(properties, c(foundation_props, water_ret_props))

  if (length(chem_props) > 0) {
    if (verbose) log_message("INFO", "Phase 3: Processing chemical properties", category = "Infilling")

    for (prop in chem_props) {
      if (verbose) log_message("DEBUG", paste("Processing", prop), category = "Infilling")

      result <- process_single_property(df, prop, max_depth, verbose)
      df <- result$data
      processing_log[[prop]] <- result$log
    }
  }

  # Step 4: Apply filtering if requested
  if (remove_unsuitable || remove_incomplete) {
    df <- apply_data_filtering(df, remove_unsuitable, remove_incomplete,
                               required_properties, properties, verbose)
  }

  # Step 5: Generate comprehensive summary
  end_time <- Sys.time()

  if (verbose) {
    generate_processing_summary(processing_log, df, start_time, end_time)
  }

  # Add processing metadata
  attr(df, "processing_log") <- processing_log
  attr(df, "processing_time") <- as.numeric(difftime(end_time, start_time, units = "secs"))
  attr(df, "properties_processed") <- properties

  return(df)
}

#' Main Soil Property Infilling Function
#'
#' Enhanced version of the core infilling function with comprehensive recovery strategies
#' and automatic exclusion of unsuitable horizons (R, Cr, O horizons).
#'
#' @param df Input soil data frame
#' @param property_name Name of the property to infill
#' @param property_config Optional property configuration
#' @param max_depth Maximum depth for infilling (default: 250 cm)
#' @param verbose Whether to provide detailed progress messages
#'
#' @return Data frame with infilled property data
#'
#' @export
infill_soil_property <- function(df,
                                 property_name,
                                 property_config = NULL,
                                 max_depth = 250,
                                 verbose = getOption("ssurgo.verbose", FALSE)) {

  if (!is.character(property_name) || length(property_name) != 1) {
    stop("property_name must be a single character string")
  }

  if (verbose) {
    log_message("INFO", paste("Starting infilling for property:", property_name), category = "Infilling")
  }

  # Special handling for RFV
  if (property_name == "rfv") {
    return(infill_rfv_property_integrated(df, max_depth, verbose))
  }

  # Validate property column exists
  r_col <- paste0(property_name, "_r")
  if (!r_col %in% names(df)) {
    stop(paste("Column", r_col, "not found in data frame"))
  }

  # Initialize essential columns
  df <- ensure_infilling_columns(df, verbose = FALSE)

  # Clean and validate property data
  if (verbose) log_message("DEBUG", "Cleaning property data", category = "Infilling")

  clean_result <- clean_property_data(df, property_name, verbose = FALSE)
  df <- clean_result$data

  # Detect unsuitable horizons using Module 0
  df$unsuitable_horizon <- is_unsuitable(df, hzname_col = "hzname")

  # Get property configuration
  if (is.null(property_config)) {
    property_config <- get_default_property_config(property_name)
  }

  # Apply depth constraints
  if ("hzdepb_r" %in% names(df)) {
    depth_mask <- df$hzdepb_r <= max_depth
  } else {
    depth_mask <- rep(TRUE, nrow(df))
  }

  # Identify problematic horizons (missing data in suitable horizons within depth limit)
  problematic_mask <- is.na(df[[r_col]]) &
    !df$unsuitable_horizon &
    depth_mask

  if (!any(problematic_mask)) {
    if (verbose) {
      log_message("INFO", "No missing values found in suitable horizons", category = "Infilling")
    }

    # Still apply range infilling
    df <- infill_property_range_values(df, property_name, property_config)
    return(df)
  }

  # Process by component groups if available
  group_col <- determine_grouping_column(df)

  if (!is.null(group_col)) {
    if (verbose) {
      log_message("DEBUG", paste("Processing by groups using column:", group_col), category = "Infilling")
    }

    df <- df |>
      dplyr::group_by(!!rlang::sym(group_col)) |>
      dplyr::group_modify(~ process_property_group(.x, property_name, problematic_mask,
                                                   property_config, max_depth, verbose)) |>
      dplyr::ungroup()
  } else {
    # Process entire dataset as one group
    df <- process_property_group(df, property_name, problematic_mask,
                                 property_config, max_depth, verbose)
  }

  # Strategy 4: Cross-component interpolation - needs to see data from OTHER
  # groups (e.g. other cokeys), so it runs at the whole-dataset level, not
  # inside the per-group pass above.
  if (verbose) {
    log_message("DEBUG", "Applying cross-component interpolation", category = "Infilling")
  }
  df <- cross_component_property_interpolation(df, r_col)

  # Strategy 5: Related-property estimation (row-independent pedological
  # relationships) - also applied at the whole-dataset level for consistency.
  if (verbose) {
    log_message("DEBUG", "Applying related-property estimation", category = "Infilling")
  }
  df <- related_property_estimation(df, property_name, property_config)

  # Strategy 6: Fallback to per-group statistics for anything still missing,
  # applied last so Strategies 4/5 get first chance at genuinely recoverable cells.
  if (!is.null(group_col)) {
    df <- df |>
      dplyr::group_by(!!rlang::sym(group_col)) |>
      dplyr::group_modify(~ process_property_group_fallback(.x, property_name,
                                                            property_config, max_depth)) |>
      dplyr::ungroup()
  } else {
    df <- process_property_group_fallback(df, property_name, property_config, max_depth)
  }

  # Final range infilling for _l and _h values
  df <- infill_property_range_values(df, property_name, property_config)

  if (verbose) {
    final_missing <- sum(is.na(df[[r_col]]) & !df$unsuitable_horizon & depth_mask)
    log_message("INFO", paste("Infilling complete. Remaining missing values:", final_missing), category = "Infilling")
  }

  return(df)
}

# ==============================================================================
# 2. HORIZON SUITABILITY AND DATA CLEANING
# ==============================================================================

#' Enhanced Property Data Cleaning
#'
#' Advanced data cleaning with statistical outlier detection, string parsing,
#' and quality reporting optimized for soil property infilling.
#'
#' @param df Input data frame
#' @param property_name Name of the property to clean
#' @param validation_config Optional validation configuration
#' @param generate_report Whether to generate cleaning report
#' @param verbose Whether to provide progress messages
#'
#' @return List containing cleaned data and optional report
#'
#' @export
clean_property_data <- function(df,
                                         property_name,
                                         validation_config = NULL,
                                         generate_report = FALSE,
                                         verbose = FALSE) {

  if (verbose) {
    log_message("DEBUG", paste("Cleaning property data for:", property_name), category = "DataCleaning")
  }

  start_time <- Sys.time()
  original_df <- df

  property_cols <- paste0(property_name, c("_r", "_l", "_h"))
  available_cols <- intersect(property_cols, names(df))

  if (length(available_cols) == 0) {
    if (verbose) {
      log_message("WARN", paste("No property columns found for", property_name), category = "DataCleaning")
    }
    return(list(data = df))
  }

  # Initialize tracking
  cleaning_actions <- list()
  type_conversions <- list()
  outliers_detected <- list()

  for (col in available_cols) {
    original_vals <- df[[col]]
    original_type <- class(original_vals)[1]

    # Advanced string parsing for character/factor data
    if (is.character(original_vals) || is.factor(original_vals)) {
      parsing_result <- advanced_string_parser_vectorized(original_vals)
      df[[col]] <- parsing_result$values

      if (verbose && parsing_result$metadata$success_rate < 0.9) {
        log_message("WARN", paste("Low parsing success rate for", col, ":",
                                  round(parsing_result$metadata$success_rate * 100, 1), "%"),
                    category = "DataCleaning")
      }
    } else {
      # Standard numeric conversion
      df[[col]] <- vectorized_type_conversion(original_vals)
    }

    # Track type conversions
    type_conversions[[col]] <- list(
      from = original_type,
      to = "numeric",
      success_rate = sum(is.finite(df[[col]])) / length(df[[col]])
    )

    # Handle infinite values
    infinite_mask <- !is.finite(df[[col]])
    df[[col]][infinite_mask] <- NA

    if (sum(infinite_mask) > 0) {
      cleaning_actions[[paste0(col, "_infinite")]] <- sum(infinite_mask)
    }

    # Statistical outlier detection for _r columns only
    if (col == paste0(property_name, "_r")) {
      outlier_result <- detect_statistical_outliers_soil_aware(df[[col]], df, property_name)

      if (outlier_result$statistics$total_outliers > 0) {
        df[[col]][outlier_result$outliers] <- NA
        outliers_detected[[col]] <- outlier_result
        cleaning_actions[[paste0(col, "_outliers")]] <- outlier_result$statistics$total_outliers

        if (verbose) {
          log_message("DEBUG", paste("Removed", outlier_result$statistics$total_outliers,
                                     "outliers from", col), category = "DataCleaning")
        }
      }

      # Apply basic range limits
      df[[col]] <- apply_basic_range_limits(df[[col]], property_name)
    }
  }

  # Generate report if requested
  if (generate_report) {
    end_time <- Sys.time()

    report <- list(
      property = property_name,
      timestamp = end_time,
      processing_time = difftime(end_time, start_time, units = "secs"),
      rows_processed = nrow(df),
      cleaning_actions = cleaning_actions,
      type_conversions = type_conversions,
      outliers_detected = outliers_detected
    )

    return(list(data = df, report = report))
  }

  return(list(data = df))
}

# ==============================================================================
# 3. PROPERTY CONFIGURATION AND VALIDATION
# ==============================================================================

#' Get Default Property Configuration
#'
#' Returns default configuration parameters for common soil properties
#' used in the infilling process.
#'
#' @param property_name Name of the soil property
#'
#' @return List with property configuration parameters
#'
#' @export
get_default_property_config <- function(property_name) {

  # Texture properties
  if (property_name %in% c('sandtotal', 'claytotal', 'silttotal')) {
    return(list(
      type = 'texture',
      units = '%',
      typical_range = c(0, 100),
      fallback_range = 8,
      sum_constraint = 100,
      related_properties = c('sandtotal', 'claytotal', 'silttotal')
    ))
  }

  # Bulk density
  else if (property_name == 'dbovendry') {
    return(list(
      type = 'bulk_density',
      units = 'g/cm^3',
      typical_range = c(0.8, 2.2),
      fallback_range = 0.15,
      horizon_effects = TRUE
    ))
  }

  # Water retention
  else if (property_name %in% c('wthirdbar', 'wfifteenbar')) {
    return(list(
      type = 'water_retention',
      units = '%',
      typical_range = c(0, 60),
      fallback_range = 3,
      clay_dependent = TRUE,
      related_properties = c('wthirdbar', 'wfifteenbar')
    ))
  }

  # Rock fragment volume
  else if (property_name == 'rfv') {
    return(list(
      type = 'rock_fragments',
      units = '%',
      typical_range = c(0.01, 85),
      fallback_range = 5,
      zero_handling = 'special'
    ))
  }

  # Cation Exchange Capacity
  else if (property_name %in% c('cec7', 'cec82', 'ecec')) {
    return(list(
      type = 'cec',
      units = 'cmol(+)/kg',
      typical_range = c(0, 100),
      fallback_range = 5,
      clay_dependent = TRUE,
      om_dependent = TRUE,
      horizon_effects = TRUE,
      related_properties = c('claytotal', 'om')
    ))
  }

  # pH
  else if (property_name %in% c('ph1to1h2o', 'ph01mcacl2', 'phh2o', 'phkcl')) {
    return(list(
      type = 'ph',
      units = 'pH units',
      typical_range = c(3.0, 10.0),
      fallback_range = 0.5,
      horizon_effects = TRUE,
      parent_material_effects = TRUE
    ))
  }

  # Organic Matter/Carbon
  else if (property_name %in% c('om', 'oc', 'ompc', 'soc')) {
    return(list(
      type = 'organic_matter',
      units = ifelse(property_name %in% c('om', 'ompc'), '%', 'g/kg'),
      typical_range = c(0, 50),
      fallback_range = 1.5,
      horizon_effects = TRUE,
      depth_dependent = TRUE,
      related_properties = c('cec7', 'ph1to1h2o')
    ))
  }

  # Generic property (fallback)
  else {
    return(list(
      type = 'generic',
      units = 'unknown',
      typical_range = NULL,
      fallback_range = 5,
      horizon_effects = FALSE
    ))
  }
}

#' Create Custom Property Configuration
#'
#' Creates a custom configuration for specialized soil properties.
#'
#' @param property_name Name of the property
#' @param property_type Type of property
#' @param units Units of measurement
#' @param typical_range Numeric vector c(min, max) of typical values
#' @param fallback_range Default spread for range calculation
#' @param related_properties Character vector of related property names
#' @param special_options List of additional options
#'
#' @return Property configuration list
#'
#' @export
create_custom_property_config <- function(property_name,
                                          property_type = "generic",
                                          units = "unknown",
                                          typical_range = NULL,
                                          fallback_range = 5,
                                          related_properties = NULL,
                                          special_options = NULL) {

  # Validate inputs
  if (!is.character(property_name) || length(property_name) != 1) {
    stop("property_name must be a single character string")
  }

  valid_types <- c("texture", "bulk_density", "water_retention", "rock_fragments",
                   "cec", "ph", "organic_matter", "generic")

  if (!property_type %in% valid_types) {
    warning(paste("property_type should be one of:", paste(valid_types, collapse = ", ")))
  }

  # Build configuration
  config <- list(
    type = property_type,
    units = units,
    fallback_range = fallback_range
  )

  # Add typical range if provided
  if (!is.null(typical_range)) {
    if (!is.numeric(typical_range) || length(typical_range) != 2) {
      stop("typical_range must be a numeric vector of length 2: c(min, max)")
    }
    if (typical_range[1] >= typical_range[2]) {
      stop("typical_range[1] must be less than typical_range[2]")
    }
    config$typical_range <- typical_range
  }

  # Add related properties
  if (!is.null(related_properties)) {
    if (!is.character(related_properties)) {
      stop("related_properties must be a character vector")
    }
    config$related_properties <- related_properties
  }

  # Add special options
  if (!is.null(special_options)) {
    if (!is.list(special_options)) {
      stop("special_options must be a list")
    }
    config <- c(config, special_options)
  }

  return(config)
}

#' Validate a Property Configuration
#'
#' Validates that a property configuration list has the required fields and
#' that they are well-formed.
#'
#' @param config Property configuration list
#' @param property_name Name of the property (for error messages)
#'
#' @return `TRUE` if valid; stops with an error if invalid
#'
#' @export
validate_property_config <- function(config, property_name) {

  if (!is.list(config)) {
    stop(paste("Configuration for", property_name, "must be a list"))
  }

  required_fields <- c("type", "units", "fallback_range")
  missing_fields <- setdiff(required_fields, names(config))
  if (length(missing_fields) > 0) {
    stop(paste("Configuration for", property_name, "missing required fields:",
               paste(missing_fields, collapse = ", ")))
  }

  if (!is.numeric(config$fallback_range) || config$fallback_range <= 0) {
    stop(paste("fallback_range for", property_name, "must be a positive number"))
  }

  if (!is.null(config$typical_range)) {
    if (!is.numeric(config$typical_range) || length(config$typical_range) != 2) {
      stop(paste("typical_range for", property_name, "must be a numeric vector of length 2"))
    }
    if (config$typical_range[1] >= config$typical_range[2]) {
      stop(paste("typical_range for", property_name, "must have min < max"))
    }
  }

  return(TRUE)
}

#' Categorize Rock Fragment Volume Values
#'
#' Categorizes a rock fragment volume (RFV) percentage into standard ranges.
#'
#' @param rfv_value Rock fragment volume percentage
#'
#' @return Character string indicating category: "none", "low", "moderate",
#'   "high", "very_high", or "extreme"
#'
#' @export
get_rfv_range_category <- function(rfv_value) {

  if (is.na(rfv_value) || rfv_value <= 0) {
    return("none")
  } else if (rfv_value <= 5) {
    return("low")
  } else if (rfv_value <= 15) {
    return("moderate")
  } else if (rfv_value <= 35) {
    return("high")
  } else if (rfv_value <= 60) {
    return("very_high")
  } else {
    return("extreme")
  }
}

#' Apply Property-Specific Constraints
#'
#' Applies typical-range and property-type-specific constraints to a numeric
#' vector of values.
#'
#' @param values Numeric vector of values to constrain
#' @param property_config Property configuration list
#'
#' @return Constrained values
#'
#' @export
apply_property_constraints <- function(values, property_config) {

  if (is.null(values) || all(is.na(values))) {
    return(values)
  }

  if (!is.null(property_config$typical_range)) {
    values <- pmax(values, property_config$typical_range[1], na.rm = TRUE)
    values <- pmin(values, property_config$typical_range[2], na.rm = TRUE)
  }

  if (property_config$type == "texture") {
    values <- pmax(values, 0, na.rm = TRUE)
    values <- pmin(values, 100, na.rm = TRUE)
  } else if (property_config$type == "ph") {
    values <- pmax(values, 0, na.rm = TRUE)
    values <- pmin(values, 14, na.rm = TRUE)
  } else if (property_config$type == "rock_fragments") {
    values <- pmax(values, 0.01, na.rm = TRUE)
    values <- pmin(values, 95, na.rm = TRUE)
  } else {
    values <- pmax(values, 0, na.rm = TRUE)
  }

  return(values)
}

#' Create a Validation Configuration
#'
#' Creates an empty, pluggable validation-rule configuration that
#' [add_range_rule()] and [add_relationship_rule()] can add rules to, and
#' [apply_validation_rules()] can apply.
#'
#' @return A validation configuration list with `range_rules`,
#'   `relationship_rules`, `conditional_rules`, and `custom_rules` slots
#'
#' @export
create_validation_config <- function() {
  list(
    range_rules = list(),
    relationship_rules = list(),
    conditional_rules = list(),
    custom_rules = list()
  )
}

#' Add a Range Validation Rule
#'
#' Adds a rule to a validation configuration flagging values of `property`
#' outside `[min_val, max_val]`.
#'
#' @param config Validation configuration from [create_validation_config()]
#' @param property Name of the property this rule applies to
#' @param min_val Minimum acceptable value
#' @param max_val Maximum acceptable value
#' @param severity Rule severity, e.g. `"error"` or `"warning"`
#'
#' @return Updated validation configuration
#'
#' @export
add_range_rule <- function(config, property, min_val, max_val, severity = "error") {
  config$range_rules[[length(config$range_rules) + 1]] <- list(
    property = property,
    min = min_val,
    max = max_val,
    severity = severity
  )
  return(config)
}

#' Add a Relationship Validation Rule
#'
#' Adds a cross-property relationship rule (e.g. a sum-to-100 texture
#' constraint) to a validation configuration. Note: as in the original
#' reference implementation, relationship rules are recorded on the config
#' but are not yet consumed by [apply_validation_rules()] (which currently
#' only enforces `range_rules`) - this mirrors the legacy behavior exactly,
#' not a bug introduced here.
#'
#' @param config Validation configuration from [create_validation_config()]
#' @param properties Character vector of properties involved in the relationship
#' @param relationship_type Type of relationship, e.g. `"sum"`
#' @param expected_sum Expected sum, when `relationship_type == "sum"`
#' @param tolerance Allowed tolerance around `expected_sum`
#'
#' @return Updated validation configuration
#'
#' @export
add_relationship_rule <- function(config, properties, relationship_type,
                                  expected_sum = NULL, tolerance = 0.1) {
  config$relationship_rules[[length(config$relationship_rules) + 1]] <- list(
    properties = properties,
    type = relationship_type,
    expected_sum = expected_sum,
    tolerance = tolerance
  )
  return(config)
}

#' Apply Validation Rules to Values
#'
#' Applies a validation configuration's range rules to a numeric vector,
#' flagging out-of-range values.
#'
#' @param values Numeric vector of values to validate
#' @param property_name Name of the property `values` represents
#' @param validation_config Validation configuration from [create_validation_config()]
#'
#' @return List with `violations` (logical vector, one per value) and
#'   `rules_applied` (count of range rules evaluated)
#'
#' @export
apply_validation_rules <- function(values, property_name, validation_config) {

  violations <- rep(FALSE, length(values))

  for (rule in validation_config$range_rules) {
    if (rule$property == property_name) {
      range_violations <- (values < rule$min | values > rule$max) & !is.na(values)
      violations <- violations | range_violations
    }
  }

  return(list(
    violations = violations,
    rules_applied = length(validation_config$range_rules)
  ))
}

#' Summarize Unsuitable Horizons
#'
#' Restores the reporting half of the legacy `filter_unsuitable_horizons()`
#' (its flagging half is already handled by [is_unsuitable()]): logs and
#' returns a summary of which horizon types were excluded from infilling.
#' Uses [log_message()] (this package's established logging convention)
#' rather than the legacy `cat()`, a deliberate improvement consistent with
#' how every other diagnostic message in this package is emitted.
#'
#' @param df Data frame with an `unsuitable_horizon` logical column (see
#'   [is_unsuitable()])
#' @param hzname_col Name of the horizon-name column
#'
#' @return List with `n_unsuitable` (count) and `horizon_types` (unique
#'   excluded horizon names)
#'
#' @export
summarize_unsuitable_horizons <- function(df, hzname_col = "hzname") {

  if (!"unsuitable_horizon" %in% names(df)) {
    stop("df must have an 'unsuitable_horizon' column (see is_unsuitable())")
  }

  n_unsuitable <- sum(df$unsuitable_horizon, na.rm = TRUE)
  horizon_types <- character(0)

  if (n_unsuitable > 0 && hzname_col %in% names(df)) {
    horizon_types <- unique(df[[hzname_col]][df$unsuitable_horizon])
    horizon_types <- horizon_types[!is.na(horizon_types)]

    log_message("INFO", paste("Found", n_unsuitable, "unsuitable horizons for infilling:",
                              paste(horizon_types, collapse = ", ")), category = "Infilling")
  }

  return(list(n_unsuitable = n_unsuitable, horizon_types = horizon_types))
}

# ==============================================================================
# 4. RANGE VALUE INFILLING
# ==============================================================================

#' Infill Property Range Values
#'
#' Infills missing _l and _h values using data-driven approaches and
#' pedological knowledge, excluding unsuitable horizons.
#'
#' @param df Input data frame
#' @param property_name Name of the property
#' @param property_config Property configuration list
#'
#' @return Data frame with infilled range values
#'
#' @export
infill_property_range_values <- function(df, property_name, property_config) {

  # Learn from existing complete ranges (only suitable horizons)
  learned_ranges <- learn_property_ranges(df, property_name, property_config)

  # Get contextual ranges based on pedological knowledge
  context_ranges <- get_property_contextual_ranges(df, property_name, property_config)

  r_col <- paste0(property_name, "_r")
  l_col <- paste0(property_name, "_l")
  h_col <- paste0(property_name, "_h")

  # Initialize range columns if they don't exist
  if (!l_col %in% names(df)) df[[l_col]] <- NA
  if (!h_col %in% names(df)) df[[h_col]] <- NA

  # Only infill suitable horizons
  suitable_mask <- if ("unsuitable_horizon" %in% names(df)) {
    !df$unsuitable_horizon
  } else {
    rep(TRUE, nrow(df))
  }

  # PERF: calculate_property_lower_bound()/calculate_property_upper_bound() (and the
  # get_contextual_spread() they call) only ever read 3 columns off `row` - r_col, "hzname", and
  # "hzdepb_r" - but apply(df[mask, ], 1, ...) forced as.matrix() to coerce EVERY column of df
  # (typically dozens) to a single character matrix first, one row per horizon needing infilling
  # (see PERFORMANCE_IMPROVEMENT_PLAN.md Tier 2). Pre-extracting just those 3 columns as plain
  # vectors, then building a small named-list "row" per element (list `[[`/`names()` behave
  # identically to a data.frame row for these helpers' purposes) avoids both the wasted columns
  # and the coercion - and is faster than slicing single-row data.frame subsets per element too
  # (data.frame `[.data.frame` has enough per-call overhead of its own to erode the savings from
  # dropping unused columns; confirmed empirically - see PERFORMANCE_IMPROVEMENT_PLAN.md's
  # benchmark note for this item). This also fixes a real latent bug the old coercion caused: with
  # `row` coerced to all-character, calculate_property_lower_bound()/upper_bound()'s
  # `if (depth <= 30)` depth-zone check (depth = row[['hzdepb_r']], never wrapped in
  # as.numeric()) was comparing STRINGS (e.g. "5" <= "30" is FALSE lexicographically, even though
  # 5 <= 30 numerically) - silently misclassifying depth zones for any horizon whose numeric
  # hzdepb_r string didn't happen to sort the same as its numeric value.
  r_vec <- df[[r_col]]
  hzname_vec <- if ("hzname" %in% names(df)) df$hzname else NULL
  depth_vec <- if ("hzdepb_r" %in% names(df)) df$hzdepb_r else NULL

  make_row <- function(i) {
    row <- stats::setNames(list(r_vec[i]), r_col)
    if (!is.null(hzname_vec)) row[["hzname"]] <- hzname_vec[i]
    if (!is.null(depth_vec)) row[["hzdepb_r"]] <- depth_vec[i]
    row
  }

  # Infill missing _l values
  missing_l_mask <- is.na(df[[l_col]]) & !is.na(df[[r_col]]) & suitable_mask

  if (any(missing_l_mask)) {
    df[missing_l_mask, l_col] <- vapply(which(missing_l_mask), function(i) {
      calculate_property_lower_bound(make_row(i), property_name, learned_ranges,
                                     context_ranges, property_config)
    }, numeric(1))
  }

  # Infill missing _h values
  missing_h_mask <- is.na(df[[h_col]]) & !is.na(df[[r_col]]) & suitable_mask

  if (any(missing_h_mask)) {
    df[missing_h_mask, h_col] <- vapply(which(missing_h_mask), function(i) {
      calculate_property_upper_bound(make_row(i), property_name, learned_ranges,
                                     context_ranges, property_config)
    }, numeric(1))
  }

  # Apply bounds checking
  if (!is.null(property_config$typical_range)) {
    df[[l_col]][suitable_mask] <- pmax(df[[l_col]][suitable_mask],
                                       property_config$typical_range[1], na.rm = TRUE)
    df[[h_col]][suitable_mask] <- pmin(df[[h_col]][suitable_mask],
                                       property_config$typical_range[2], na.rm = TRUE)
  }

  # Ensure logical order: _l <= _r <= _h
  df[[l_col]] <- as.numeric(df[[l_col]])
  df[[r_col]] <- as.numeric(df[[r_col]])
  df[[h_col]] <- as.numeric(df[[h_col]])

  df[[l_col]][suitable_mask] <- pmin(df[[l_col]][suitable_mask],
                                     df[[r_col]][suitable_mask], na.rm = TRUE)
  df[[h_col]][suitable_mask] <- pmax(df[[h_col]][suitable_mask],
                                     df[[r_col]][suitable_mask], na.rm = TRUE)

  return(df)
}

#' Learn Property Ranges from Data
#'
#' Learns typical ranges from existing complete data for any property,
#' only using suitable horizons for learning.
#'
#' @param df Input data frame
#' @param property_name Name of the property
#' @param property_config Property configuration
#'
#' @return List of learned ranges by context
#'
#' @export
learn_property_ranges <- function(df, property_name, property_config) {

  r_col <- paste0(property_name, "_r")
  l_col <- paste0(property_name, "_l")
  h_col <- paste0(property_name, "_h")

  # Ensure columns exist and are numeric
  for (col in c(r_col, l_col, h_col)) {
    if (col %in% names(df)) {
      df[[col]] <- as.numeric(df[[col]])
    }
  }

  # Get complete records from suitable horizons only
  complete_mask <- !is.na(df[[r_col]]) & !is.na(df[[l_col]]) & !is.na(df[[h_col]]) &
    is.finite(df[[r_col]]) & is.finite(df[[l_col]]) & is.finite(df[[h_col]])

  # Exclude unsuitable horizons
  if ("unsuitable_horizon" %in% names(df)) {
    suitable_mask <- !df$unsuitable_horizon
    complete_mask <- complete_mask & suitable_mask
  }

  complete_data <- df[complete_mask, ]

  if (nrow(complete_data) == 0) {
    fallback_range <- ifelse(is.null(property_config$fallback_range), 5, property_config$fallback_range)
    return(list(default_spread = as.numeric(fallback_range)))
  }

  # Calculate actual spreads
  lower_spreads <- complete_data[[r_col]] - complete_data[[l_col]]
  upper_spreads <- complete_data[[h_col]] - complete_data[[r_col]]

  ranges_by_context <- list()

  # Learn ranges by horizon type
  if ('hzname' %in% names(complete_data)) {
    hznames <- unique(complete_data$hzname[!is.na(complete_data$hzname)])

    for (hzname in hznames) {
      hz_data <- complete_data[complete_data$hzname == hzname & !is.na(complete_data$hzname), ]

      if (nrow(hz_data) >= 3) {
        hz_lower <- hz_data[[r_col]] - hz_data[[l_col]]
        hz_upper <- hz_data[[h_col]] - hz_data[[r_col]]

        ranges_by_context[[paste0('hzname_', hzname)]] <- list(
          lower_spread_median = median(hz_lower, na.rm = TRUE),
          upper_spread_median = median(hz_upper, na.rm = TRUE),
          sample_size = nrow(hz_data)
        )
      }
    }
  }

  # Learn ranges by depth zone
  if ('hzdepb_r' %in% names(complete_data)) {
    depth_zones <- list(
      surface = c(0, 30),
      subsurface = c(30, 100),
      deep = c(100, 200)
    )

    for (zone_name in names(depth_zones)) {
      zone_range <- depth_zones[[zone_name]]
      zone_data <- complete_data[complete_data$hzdepb_r > zone_range[1] &
                                   complete_data$hzdepb_r <= zone_range[2], ]

      if (nrow(zone_data) >= 3) {
        zone_lower <- zone_data[[r_col]] - zone_data[[l_col]]
        zone_upper <- zone_data[[h_col]] - zone_data[[r_col]]

        ranges_by_context[[paste0('depth_', zone_name)]] <- list(
          lower_spread_median = median(zone_lower, na.rm = TRUE),
          upper_spread_median = median(zone_upper, na.rm = TRUE),
          sample_size = nrow(zone_data)
        )
      }
    }
  }

  # Overall dataset statistics
  ranges_by_context[['overall']] <- list(
    lower_spread_median = median(lower_spreads, na.rm = TRUE),
    upper_spread_median = median(upper_spreads, na.rm = TRUE),
    sample_size = nrow(complete_data)
  )

  return(ranges_by_context)
}

#' Get Property Contextual Ranges
#'
#' Returns property-specific contextual ranges based on soil science knowledge.
#'
#' @param df Input data frame
#' @param property_name Name of the property
#' @param property_config Property configuration
#'
#' @return List of contextual ranges
#'
#' @export
get_property_contextual_ranges <- function(df, property_name, property_config) {

  contextual_ranges <- list()

  # Property-specific knowledge
  if (property_config$type == 'texture') {
    contextual_ranges[[property_name]] <- switch(property_name,
                                                 sandtotal = list(typical_spread = 12, min_spread = 5, max_spread = 25),
                                                 claytotal = list(typical_spread = 8, min_spread = 3, max_spread = 20),
                                                 silttotal = list(typical_spread = 10, min_spread = 4, max_spread = 22),
                                                 list(typical_spread = 10, min_spread = 5, max_spread = 20)
    )
  }

  else if (property_config$type == 'bulk_density') {
    contextual_ranges <- list(
      A_horizons = list(typical_spread = 0.15, range = c(0.8, 1.6)),
      B_horizons = list(typical_spread = 0.12, range = c(1.0, 1.8)),
      C_horizons = list(typical_spread = 0.10, range = c(1.2, 2.0))
    )
  }

  else if (property_config$type == 'water_retention') {
    contextual_ranges[[property_name]] <- switch(property_name,
                                                 wthirdbar = list(typical_spread = 3, clay_factor = 0.4),
                                                 wfifteenbar = list(typical_spread = 2, clay_factor = 0.3),
                                                 list(typical_spread = 2.5, clay_factor = 0.35)
    )
  }

  else if (property_config$type == 'rock_fragments') {
    contextual_ranges[['rfv_ranges']] <- list(
      low = list(range = c(0.01, 5), typical_spread = 2),
      moderate = list(range = c(5, 15), typical_spread = 4),
      high = list(range = c(15, 35), typical_spread = 6),
      very_high = list(range = c(35, 60), typical_spread = 10),
      extreme = list(range = c(60, 85), typical_spread = 12)
    )
  }

  else {
    # Generic property
    contextual_ranges[['generic']] <- list(
      typical_spread = property_config$fallback_range,
      percentage_based = TRUE
    )
  }

  return(contextual_ranges)
}

# ==============================================================================
# 5. PROPERTY DATA RECOVERY STRATEGIES
# ==============================================================================

#' Infill Missing Property Data
#'
#' Applies the "within-group" recovery strategies, in order of reliability,
#' only using suitable horizons: horizon-name matching, depth-weighted
#' averaging, and within-component interpolation. This function is called
#' once per grouping unit (typically per `cokey`) via [process_property_group()],
#' so it cannot see data from other groups - the group-spanning strategies
#' (cross-component interpolation, related-property estimation, and the final
#' group-mean fallback) are applied afterward at the whole-dataset level by
#' [infill_soil_property()] (see [cross_component_property_interpolation()],
#' [related_property_estimation()], and `apply_group_fallback_mean()`).
#'
#' @param group Data frame group to process
#' @param property_name Name of the property to infill
#' @param problematic_mask Logical vector for problematic horizons
#' @param property_config Property configuration
#'
#' @return Data frame with infilled property data
#'
#' @export
infill_missing_property_data <- function(group, property_name, problematic_mask, property_config) {

  property_col <- paste0(property_name, "_r")
  has_depth <- all(c('hzdept_r', 'hzdepb_r') %in% names(group))

  # Validate problematic_mask length matches group
  if (length(problematic_mask) != nrow(group)) {
    warning(paste("Problematic mask length mismatch. Recalculating for", property_name))

    # Recalculate if mismatch
    suitable_mask <- if ("unsuitable_horizon" %in% names(group)) {
      !group$unsuitable_horizon
    } else {
      rep(TRUE, nrow(group))
    }

    problematic_mask <- is.na(group[[property_col]]) & suitable_mask
  }

  # Strategy 1: Horizon name matching (highest priority)
  if ('hzname' %in% names(group)) {
    group <- horizon_name_property_infill(group, property_col, problematic_mask)
  }

  # Check remaining missing data
  still_missing <- is.na(group[[property_col]]) & problematic_mask
  if (!any(still_missing)) return(group)

  # Strategy 2: Depth-weighted averaging
  if (has_depth) {
    group <- depth_weighted_property_infill(group, property_col, problematic_mask)
  }

  still_missing <- is.na(group[[property_col]]) & problematic_mask
  if (!any(still_missing)) return(group)

  # Strategy 3: Within-component interpolation
  if (has_depth) {
    group <- within_component_property_interpolation(group, property_col)
  }

  return(group)
}

#' Apply Group-Mean Fallback
#'
#' Strategy 6 (last resort) of the property recovery hierarchy: fills any
#' still-missing suitable-horizon values with a depth-weighted (or plain)
#' mean computed from the suitable horizons of `group`. Extracted as its own
#' function so it can be applied after the group-spanning Strategy 4/5 passes
#' (cross-component interpolation, related-property estimation) have already
#' had a chance to fill values - see [infill_soil_property()].
#'
#' @param group Data frame group to process
#' @param property_col Name of the property column to infill
#' @param problematic_mask Logical vector for problematic horizons
#'
#' @return Data frame with fallback-filled values
apply_group_fallback_mean <- function(group, property_col, problematic_mask) {

  has_depth <- all(c('hzdept_r', 'hzdepb_r') %in% names(group))
  remaining_missing <- is.na(group[[property_col]]) & problematic_mask

  if (any(remaining_missing)) {
    suitable_mask <- if ("unsuitable_horizon" %in% names(group)) {
      !group$unsuitable_horizon
    } else {
      rep(TRUE, nrow(group))
    }

    suitable_group <- group[suitable_mask, ]

    if (has_depth && nrow(suitable_group) > 0) {
      depth_weighted_mean <- calculate_depth_weighted_mean(suitable_group, property_col)
    } else if (nrow(suitable_group) > 0) {
      depth_weighted_mean <- mean(suitable_group[[property_col]], na.rm = TRUE)
    } else {
      depth_weighted_mean <- NA
    }

    if (!is.na(depth_weighted_mean)) {
      group[[property_col]][remaining_missing] <- depth_weighted_mean

      # Update tracking
      if (!"infill_method" %in% names(group)) {
        group$infill_method <- ""
      }

      group$infill_method[remaining_missing] <- paste0(group$infill_method[remaining_missing],
                                                       property_col, ":group_fallback; ")
    }
  }

  return(group)
}

#' Horizon Name Property Infill
#'
#' Infills missing values using horizons with matching or similar names,
#' only using suitable horizons as sources.
#'
#' @param group Data frame group
#' @param property_col Property column name
#' @param problematic_mask Logical mask for problematic horizons
#'
#' @return Data frame with infilled values
#'
#' @export
horizon_name_property_infill <- function(group, property_col, problematic_mask) {

  if (!property_col %in% names(group) || !"hzname" %in% names(group)) {
    return(group)
  }

  # Initialize tracking
  if (!"infill_method" %in% names(group)) {
    group$infill_method <- ""
  }

  # Validate mask length
  if (length(problematic_mask) != nrow(group)) {
    warning("Mask length mismatch in horizon_name_property_infill")
    return(group)
  }

  # Determine suitable horizons for sources
  suitable_mask <- if ("unsuitable_horizon" %in% names(group)) {
    !group$unsuitable_horizon
  } else {
    rep(TRUE, nrow(group))
  }

  # Get indices needing infilling
  problematic_indices <- which(problematic_mask &
                                 is.na(group[[property_col]]) &
                                 !is.na(group$hzname))

  if (length(problematic_indices) == 0) return(group)

  # Pre-process horizon names
  standardized_hznames <- sapply(group$hzname, function(x) {
    if (is.na(x)) return(NA_character_)
    standardize_horizon_name(x)
  })

  # Find source horizons
  source_mask <- suitable_mask &
    !is.na(group[[property_col]]) &
    !is.na(standardized_hznames)

  source_indices <- which(source_mask)

  if (length(source_indices) == 0) return(group)

  # Process each problematic horizon
  for (idx in problematic_indices) {
    target_hzname <- standardized_hznames[idx]

    if (is.na(target_hzname) || target_hzname == "") next

    # Find matching horizons
    matching_values <- numeric(0)
    similarity_scores <- numeric(0)

    for (source_idx in source_indices) {
      if (source_idx == idx) next

      source_hzname <- standardized_hznames[source_idx]
      if (is.na(source_hzname) || source_hzname == "") next

      similarity <- calculate_horizon_similarity(target_hzname, source_hzname)

      if (similarity > 0.5) {
        matching_values <- c(matching_values, group[[property_col]][source_idx])
        similarity_scores <- c(similarity_scores, similarity)
      }
    }

    # Infill if matches found
    if (length(matching_values) > 0) {
      if (length(matching_values) == 1) {
        infilled_value <- matching_values[1]
      } else {
        if (all(is.finite(similarity_scores)) && all(similarity_scores > 0)) {
          infilled_value <- weighted.mean(matching_values, similarity_scores)
        } else {
          infilled_value <- mean(matching_values)
        }
      }

      group[[property_col]][idx] <- infilled_value

      # Update tracking
      method_info <- paste0(property_col, ":hzname(", target_hzname, ",n=",
                            length(matching_values), "); ")

      if (is.na(group$infill_method[idx]) || group$infill_method[idx] == "") {
        group$infill_method[idx] <- method_info
      } else {
        group$infill_method[idx] <- paste0(group$infill_method[idx], method_info)
      }
    }
  }

  return(group)
}

# ==============================================================================
# 6. SPECIAL PROPERTY HANDLING
# ==============================================================================

#' Enhanced RFV Property Infilling
#'
#' Specialized infilling for rock fragment volume with context-aware estimation
#' and integration with the main workflow.
#'
#' @param df Input data frame
#' @param max_depth Maximum depth for processing
#' @param verbose Whether to provide progress messages
#'
#' @return Data frame with infilled RFV values
#'
#' @export
infill_rfv_property_integrated <- function(df, max_depth = 250, verbose = FALSE) {

  if (verbose) {
    log_message("INFO", "Processing RFV with specialized handling", category = "RFV")
  }

  # Ensure essential columns
  df <- ensure_infilling_columns(df, verbose = FALSE)

  # Determine processing mask
  suitable_mask <- if ("unsuitable_horizon" %in% names(df)) {
    !df$unsuitable_horizon
  } else {
    rep(TRUE, nrow(df))
  }

  if ("hzdepb_r" %in% names(df)) {
    depth_mask <- df$hzdepb_r <= max_depth
    process_mask <- suitable_mask & depth_mask
  } else {
    process_mask <- suitable_mask
  }

  process_indices <- which(process_mask)

  if (length(process_indices) == 0) {
    if (verbose) {
      log_message("WARN", "No suitable horizons found for RFV processing", category = "RFV")
    }
    return(df)
  }

  if (verbose) {
    log_message("INFO", paste("Processing", length(process_indices),
                              "suitable RFV horizons"), category = "RFV")
  }

  # Process each suitable horizon
  for (i in process_indices) {
    tryCatch({
      row_list <- as.list(df[i, ])
      result_row <- impute_rfv_values(row_list)
      df[i, names(result_row)] <- result_row

      # Update tracking
      if (is.na(df$infill_method[i]) || df$infill_method[i] == "") {
        df$infill_method[i] <- "rfv_specialized; "
      } else {
        df$infill_method[i] <- paste0(df$infill_method[i], "rfv_specialized; ")
      }

    }, error = function(e) {
      if (verbose) {
        log_message("WARN", paste("RFV processing error for row", i, ":", e$message), category = "RFV")
      }
    })
  }

  # Apply range infilling
  rfv_config <- get_default_property_config("rfv")
  df <- infill_property_range_values(df, "rfv", rfv_config)

  return(df)
}

#' Water Retention Estimation using Saxton-Rawls
#'
#' Estimates missing water retention values using Saxton-Rawls equations
#' with rock fragment correction, integrated with the main workflow.
#'
#' @param df Input data frame
#' @param max_depth Maximum depth for estimation
#' @param add_ranges Whether to add range estimates
#' @param overwrite Whether to overwrite existing values
#' @param verbose Whether to provide progress messages
#'
#' @return Data frame with estimated water retention values
#'
#' @export
infill_water_retention_saxton_rawls_integrated <- function(df,
                                                           max_depth = 250,
                                                           add_ranges = TRUE,
                                                           overwrite = FALSE,
                                                           verbose = FALSE) {

  # Validate required columns
  required_cols <- c("claytotal_r", "sandtotal_r", "silttotal_r", "dbovendry_r")
  missing_cols <- setdiff(required_cols, names(df))

  if (length(missing_cols) > 0) {
    if (verbose) {
      log_message("WARN", paste("Missing required columns for water retention:",
                                paste(missing_cols, collapse = ", ")), category = "WaterRetention")
    }
    return(df)
  }

  if (verbose) {
    log_message("INFO", "Estimating water retention using Saxton-Rawls equations", category = "WaterRetention")
  }

  # Initialize columns
  if (!"wthirdbar_r" %in% names(df)) df$wthirdbar_r <- NA_real_
  if (!"wfifteenbar_r" %in% names(df)) df$wfifteenbar_r <- NA_real_

  if (add_ranges) {
    if (!"wthirdbar_l" %in% names(df)) df$wthirdbar_l <- NA_real_
    if (!"wthirdbar_h" %in% names(df)) df$wthirdbar_h <- NA_real_
    if (!"wfifteenbar_l" %in% names(df)) df$wfifteenbar_l <- NA_real_
    if (!"wfifteenbar_h" %in% names(df)) df$wfifteenbar_h <- NA_real_
  }

  # Ensure unsuitable horizon detection
  if (!"unsuitable_horizon" %in% names(df)) {
    df$unsuitable_horizon <- is_unsuitable(df)
  }

  # Determine processing mask
  suitable_mask <- !df$unsuitable_horizon

  if ("hzdepb_r" %in% names(df)) {
    depth_mask <- df$hzdepb_r <= max_depth
    process_mask <- suitable_mask & depth_mask
  } else {
    process_mask <- suitable_mask
  }

  # Check for complete input data
  complete_inputs_mask <- process_mask &
    !is.na(df$claytotal_r) & !is.na(df$sandtotal_r) &
    !is.na(df$silttotal_r) & !is.na(df$dbovendry_r) &
    is.finite(df$claytotal_r) & is.finite(df$sandtotal_r) &
    is.finite(df$silttotal_r) & is.finite(df$dbovendry_r)

  # Determine estimation needs
  if (overwrite) {
    need_fc_estimate <- complete_inputs_mask
    need_wp_estimate <- complete_inputs_mask
  } else {
    need_fc_estimate <- complete_inputs_mask & is.na(df$wthirdbar_r)
    need_wp_estimate <- complete_inputs_mask & is.na(df$wfifteenbar_r)
  }

  n_fc <- sum(need_fc_estimate)
  n_wp <- sum(need_wp_estimate)

  if (n_fc == 0 && n_wp == 0) {
    if (verbose) {
      log_message("INFO", "No water retention values need estimation", category = "WaterRetention")
    }
    return(df)
  }

  if (verbose) {
    log_message("INFO", paste("Estimating field capacity:", n_fc, "values"), category = "WaterRetention")
    log_message("INFO", paste("Estimating wilting point:", n_wp, "values"), category = "WaterRetention")
  }

  # Get optional data
  rfv_values <- if ("rfv_r" %in% names(df)) {
    ifelse(is.na(df$rfv_r) | !is.finite(df$rfv_r), 0, df$rfv_r)
  } else {
    rep(0, nrow(df))
  }

  om_values <- if ("om_r" %in% names(df)) {
    ifelse(is.na(df$om_r) | !is.finite(df$om_r), 2, df$om_r)
  } else {
    rep(2, nrow(df))
  }

  # Process each horizon needing estimation
  indices_to_process <- which(need_fc_estimate | need_wp_estimate)

  for (i in indices_to_process) {
    tryCatch({
      result <- calculate_saxton_rawls_single(
        sand_pct = df$sandtotal_r[i],
        clay_pct = df$claytotal_r[i],
        silt_pct = df$silttotal_r[i],
        bulk_density = df$dbovendry_r[i],
        rfv_pct = rfv_values[i],
        om_pct = om_values[i]
      )

      # Update values
      if (need_fc_estimate[i]) {
        df$wthirdbar_r[i] <- result$field_capacity
        if (add_ranges) {
          df$wthirdbar_l[i] <- result$field_capacity_l
          df$wthirdbar_h[i] <- result$field_capacity_h
        }
      }

      if (need_wp_estimate[i]) {
        df$wfifteenbar_r[i] <- result$wilting_point
        if (add_ranges) {
          df$wfifteenbar_l[i] <- result$wilting_point_l
          df$wfifteenbar_h[i] <- result$wilting_point_h
        }
      }

      # Update tracking
      if (!"infill_method" %in% names(df)) {
        df$infill_method <- ""
      }

      method_components <- c()
      if (need_fc_estimate[i]) method_components <- c(method_components, "wthirdbar_saxton_rawls")
      if (need_wp_estimate[i]) method_components <- c(method_components, "wfifteenbar_saxton_rawls")

      if (rfv_values[i] > 0.1) {
        method_components <- paste0(method_components, "_rfv")
      }

      method_info <- paste0(paste(method_components, collapse = ","), "; ")

      if (is.na(df$infill_method[i]) || df$infill_method[i] == "") {
        df$infill_method[i] <- method_info
      } else {
        df$infill_method[i] <- paste0(df$infill_method[i], method_info)
      }

    }, error = function(e) {
      if (verbose) {
        log_message("WARN", paste("Water retention estimation error for row", i, ":", e$message),
                    category = "WaterRetention")
      }
    })
  }

  if (verbose) {
    final_fc <- sum(!is.na(df$wthirdbar_r[need_fc_estimate]))
    final_wp <- sum(!is.na(df$wfifteenbar_r[need_wp_estimate]))
    log_message("INFO", paste("Successfully estimated - FC:", final_fc, "WP:", final_wp), category = "WaterRetention")
  }

  return(df)
}

# ==============================================================================
# 7. SUPPORTING UTILITY FUNCTIONS
# ==============================================================================

#' Auto-detect Soil Properties
#'
#' Automatically detects available soil properties in a dataset.
#'
#' @param df Input data frame
#'
#' @return Vector of detected property names
auto_detect_soil_properties <- function(df) {

  available_properties <- character(0)

  # Physical properties
  physical_props <- c("sandtotal", "claytotal", "silttotal", "dbovendry", "rfv")
  for (prop in physical_props) {
    if (paste0(prop, "_r") %in% names(df)) {
      available_properties <- c(available_properties, prop)
    }
  }

  # Chemical properties
  chemical_props <- c("cec7", "ph1to1h2o", "om")
  for (prop in chemical_props) {
    if (paste0(prop, "_r") %in% names(df)) {
      available_properties <- c(available_properties, prop)
    }
  }

  # Water retention
  water_props <- c("wthirdbar", "wfifteenbar")
  for (prop in water_props) {
    if (paste0(prop, "_r") %in% names(df)) {
      available_properties <- c(available_properties, prop)
    }
  }

  return(available_properties)
}

#' Ensure Infilling Columns
#'
#' Ensures required columns exist for the infilling process.
#'
#' @param df Input data frame
#' @param verbose Whether to provide progress messages
#'
#' @return Data frame with required columns
ensure_infilling_columns <- function(df, verbose = FALSE) {

  # Initialize infill_method column
  if (!"infill_method" %in% names(df)) {
    df$infill_method <- ""
  }

  # Ensure cokey exists
  if (!"cokey" %in% names(df)) {
    if ("compname" %in% names(df)) {
      df$cokey <- paste0("comp_", df$compname)
    } else {
      df$cokey <- paste0("missing_", seq_len(nrow(df)))
    }

    if (verbose) {
      log_message("DEBUG", "Generated cokey column", category = "Infilling")
    }
  }

  # Ensure hzname exists
  if (!"hzname" %in% names(df)) {
    df$hzname <- "Unknown"

    if (verbose) {
      log_message("DEBUG", "Generated hzname column", category = "Infilling")
    }
  }

  # Ensure numeric depth columns
  if ("hzdept_r" %in% names(df)) {
    df$hzdept_r <- as.numeric(df$hzdept_r)
  }
  if ("hzdepb_r" %in% names(df)) {
    df$hzdepb_r <- as.numeric(df$hzdepb_r)
  }

  return(df)
}

#' Determine Grouping Column
#'
#' Determines the best column for grouping soil components.
#'
#' @param df Input data frame
#'
#' @return Column name for grouping or NULL
determine_grouping_column <- function(df) {

  if ("cokey" %in% names(df)) {
    return("cokey")
  } else if ("compname" %in% names(df)) {
    return("compname")
  } else {
    return(NULL)
  }
}

#' Process Single Property
#'
#' Processes a single property with logging and error handling.
#'
#' @param df Input data frame
#' @param property_name Property to process
#' @param max_depth Maximum depth
#' @param verbose Verbose output
#'
#' @return List with processed data and log
process_single_property <- function(df, property_name, max_depth, verbose) {

  r_col <- paste0(property_name, "_r")

  # Track initial state
  initial_missing <- if (r_col %in% names(df)) {
    sum(is.na(df[[r_col]]))
  } else {
    nrow(df)
  }

  # Process the property
  if (property_name == "rfv") {
    result_df <- infill_rfv_property_integrated(df, max_depth, verbose)
  } else {
    result_df <- infill_soil_property(df, property_name, max_depth = max_depth, verbose = verbose)
  }

  # Track final state
  final_missing <- if (r_col %in% names(result_df)) {
    sum(is.na(result_df[[r_col]]))
  } else {
    nrow(result_df)
  }

  # Create log entry
  log_entry <- list(
    property = property_name,
    initial_missing = initial_missing,
    final_missing = final_missing,
    filled = initial_missing - final_missing
  )

  return(list(data = result_df, log = log_entry))
}

#' Apply Data Filtering
#'
#' Applies optional data filtering based on user preferences.
#'
#' @param df Input data frame
#' @param remove_unsuitable Remove unsuitable horizons
#' @param remove_incomplete Remove incomplete rows
#' @param required_properties Required complete properties
#' @param all_properties All processed properties
#' @param verbose Verbose output
#'
#' @return Filtered data frame
apply_data_filtering <- function(df, remove_unsuitable, remove_incomplete,
                                 required_properties, all_properties, verbose) {

  original_rows <- nrow(df)

  # Remove unsuitable horizons
  if (remove_unsuitable && "unsuitable_horizon" %in% names(df)) {
    unsuitable_count <- sum(df$unsuitable_horizon, na.rm = TRUE)
    df <- df[!df$unsuitable_horizon %in% TRUE, ]

    if (verbose) {
      log_message("INFO", paste("Removed", unsuitable_count, "unsuitable horizons"), category = "Filtering")
    }
  }

  # Remove incomplete rows
  if (remove_incomplete) {
    if (is.null(required_properties)) {
      required_properties <- all_properties
    }

    required_cols <- paste0(required_properties, "_r")
    existing_cols <- intersect(required_cols, names(df))

    if (length(existing_cols) > 0) {
      complete_mask <- complete.cases(df[, existing_cols, drop = FALSE])
      rows_before <- nrow(df)
      df <- df[complete_mask, ]

      if (verbose) {
        log_message("INFO", paste("Removed", rows_before - nrow(df), "incomplete rows"), category = "Filtering")
      }
    }
  }

  final_rows <- nrow(df)

  if (verbose) {
    log_message("INFO", paste("Filtering summary: ", original_rows, "->", final_rows,
                              "(", round(final_rows/original_rows*100, 1), "% retained)"), category = "Filtering")
  }

  return(df)
}

#' Generate Processing Summary
#'
#' Generates a comprehensive summary of the processing results.
#'
#' @param processing_log Processing log
#' @param df Final data frame
#' @param start_time Start time
#' @param end_time End time
generate_processing_summary <- function(processing_log, df, start_time, end_time) {

  log_message("INFO", "=== PROCESSING COMPLETE ===", category = "Summary")
  log_message("INFO", paste("Duration:", round(as.numeric(difftime(end_time, start_time, units = "secs")), 2), "seconds"), category = "Summary")

  if ("unsuitable_horizon" %in% names(df)) {
    unsuitable_count <- sum(df$unsuitable_horizon, na.rm = TRUE)
    log_message("INFO", paste("Unsuitable horizons excluded:", unsuitable_count), category = "Summary")
  }

  for (prop in names(processing_log)) {
    log_entry <- processing_log[[prop]]
    if (is.list(log_entry) && "filled" %in% names(log_entry)) {
      log_message("INFO", paste(sprintf("%-12s: %3d missing -> %3d missing (%3d filled)",
                                        prop, log_entry$initial_missing,
                                        log_entry$final_missing, log_entry$filled)), category = "Summary")
    }
  }
}

# ==============================================================================
# SUPPORTING FUNCTIONS FROM ORIGINAL SCRIPT
# ==============================================================================

#' Advanced String Parser (Vectorized)
#'
#' Parses character/factor soil data with advanced pattern recognition.
#'
#' @param value_strings Character vector to parse
#'
#' @return List with parsed values and metadata
advanced_string_parser_vectorized <- function(value_strings) {
  if (is.factor(value_strings)) {
    value_strings <- as.character(value_strings)
  }

  parsed_values <- rep(NA_real_, length(value_strings))
  parsing_metadata <- rep("unparseable", length(value_strings))
  confidence_scores <- rep(0, length(value_strings))

  clean_strings <- trimws(toupper(value_strings))

  for (i in seq_along(clean_strings)) {
    if (is.na(clean_strings[i]) || clean_strings[i] == "") {
      next
    }

    result <- parse_single_string_advanced(clean_strings[i])
    parsed_values[i] <- result$value
    parsing_metadata[i] <- result$type
    confidence_scores[i] <- result$confidence
  }

  return(list(
    values = parsed_values,
    metadata = list(
      parsing_types = parsing_metadata,
      confidence_scores = confidence_scores,
      success_rate = sum(!is.na(parsed_values)) / length(parsed_values)
    )
  ))
}

#' Parse Single String Advanced
#'
#' Advanced pattern recognition for individual soil property values.
#'
#' @param value_string Single string to parse
#'
#' @return List with value, type, and confidence
parse_single_string_advanced <- function(value_string) {

  # Range values: "15-20", "10 TO 15", "5-20" with an en/em dash instead of a hyphen.
  # Non-ASCII dash/comparison characters below are expressed as \uXXXX
  # escapes (not literal characters) to keep this source file ASCII-only
  # per R CMD check, with identical matching behavior at runtime.
  range_pattern <- "\\d+\\.?\\d*\\s*[-\u2013\u2014TO]\\s*\\d+\\.?\\d*"
  if (grepl(range_pattern, value_string)) {
    numbers <- as.numeric(unlist(regmatches(value_string, gregexpr("\\d+\\.?\\d*", value_string))))
    if (length(numbers) == 2 && all(is.finite(numbers))) {
      return(list(value = mean(numbers), type = "range", confidence = 0.8))
    }
  }

  # Comparison operators: "<5", ">50", using <= as a unicode less-than-or-equal sign, ">=25"
  comparison_pattern <- "^[<>\u2264\u2265]\\s*=?\\s*\\d+\\.?\\d*"
  if (grepl(comparison_pattern, value_string)) {
    number <- as.numeric(gsub("[^0-9.]", "", value_string))
    if (is.finite(number)) {
      operator <- gsub("\\d+.*", "", gsub("\\s", "", value_string))
      if (operator %in% c("<", "\u2264")) {
        return(list(value = number * 0.5, type = "upper_bound", confidence = 0.6))
      } else if (operator %in% c(">", "\u2265", ">=")) {
        return(list(value = number * 1.5, type = "lower_bound", confidence = 0.6))
      }
    }
  }

  # Qualitative terms with soil science domain knowledge
  qualitative_mapping <- list(
    "TRACE" = 0.1, "TRACES" = 0.1, "TR" = 0.1,
    "NONE" = 0.01, "ABSENT" = 0.01, "NIL" = 0.01,
    "LOW" = 15, "MODERATE" = 35, "MOD" = 35, "HIGH" = 65,
    "VERY LOW" = 5, "V LOW" = 5, "VL" = 5,
    "VERY HIGH" = 85, "V HIGH" = 85, "VH" = 85,
    "SLIGHT" = 10, "STRONG" = 60, "WEAK" = 20,
    "SANDY" = 75, "CLAYEY" = 45, "SILTY" = 50,
    "COARSE" = 70, "FINE" = 25, "MEDIUM" = 45,
    "ND" = 0.01, "BDL" = 0.01, "LOD" = 0.01,
    "SCAT" = 5, "FEW" = 8, "COMMON" = 15, "MANY" = 25
  )

  if (value_string %in% names(qualitative_mapping)) {
    return(list(value = qualitative_mapping[[value_string]], type = "qualitative", confidence = 0.4))
  }

  # Percentages: "25.5%", "30 PCT", "45 PERCENT"
  if (grepl("%|PCT|PERCENT", value_string)) {
    number <- as.numeric(gsub("[^0-9.]", "", value_string))
    if (is.finite(number)) {
      return(list(value = number, type = "percentage", confidence = 0.9))
    }
  }

  # Standard numeric extraction
  numeric_pattern <- "\\d+\\.?\\d*"
  if (grepl(numeric_pattern, value_string)) {
    number <- as.numeric(regmatches(value_string, regexpr(numeric_pattern, value_string)))
    if (is.finite(number)) {
      return(list(value = number, type = "numeric", confidence = 0.9))
    }
  }

  return(list(value = NA, type = "unparseable", confidence = 0))
}

#' Vectorized Type Conversion
#'
#' Efficient type conversion for soil property data.
#'
#' @param values Vector to convert
#'
#' @return Numeric vector
vectorized_type_conversion <- function(values) {
  if (is.factor(values)) {
    return(as.numeric(as.character(values)))
  } else if (is.character(values)) {
    return(as.numeric(gsub("[^0-9.-]", "", values)))
  } else {
    return(as.numeric(values))
  }
}

#' Statistical Outlier Detection (Soil-Aware)
#'
#' Property-specific outlier detection that respects soil science ranges.
#'
#' @param values Numeric vector to analyze
#' @param df Full dataframe for context
#' @param property_name Property being analyzed
#'
#' @return List with outlier detection results
detect_statistical_outliers_soil_aware <- function(values, df, property_name) {

  if (all(is.na(values)) || length(values) < 10) {
    return(list(
      outliers = rep(FALSE, length(values)),
      method_used = "insufficient_data",
      statistics = list(total_outliers = 0, outlier_rate = 0)
    ))
  }

  # Property-specific outlier detection
  if (property_name == "ph1to1h2o" || grepl("^ph", property_name)) {
    # pH: only flag physically impossible values
    extreme_outliers <- (values < 2.5 | values > 11.5) & !is.na(values)

    return(list(
      outliers = extreme_outliers,
      method_used = "ph_conservative",
      statistics = list(
        total_outliers = sum(extreme_outliers, na.rm = TRUE),
        outlier_rate = mean(extreme_outliers, na.rm = TRUE)
      )
    ))

  } else if (property_name %in% c("sandtotal", "claytotal", "silttotal")) {
    # Texture: only flag impossible percentages
    impossible_outliers <- (values < 0 | values > 100) & !is.na(values)

    return(list(
      outliers = impossible_outliers,
      method_used = "texture_conservative",
      statistics = list(
        total_outliers = sum(impossible_outliers, na.rm = TRUE),
        outlier_rate = mean(impossible_outliers, na.rm = TRUE)
      )
    ))

  } else if (property_name == "dbovendry") {
    # Bulk density: only flag physically impossible values
    impossible_outliers <- (values < 0.3 | values > 3.5) & !is.na(values)

    return(list(
      outliers = impossible_outliers,
      method_used = "bulk_density_conservative",
      statistics = list(
        total_outliers = sum(impossible_outliers, na.rm = TRUE),
        outlier_rate = mean(impossible_outliers, na.rm = TRUE)
      )
    ))

  } else {
    # Generic: very conservative IQR (factor = 5.0)
    Q1 <- quantile(values, 0.25, na.rm = TRUE)
    Q3 <- quantile(values, 0.75, na.rm = TRUE)
    IQR_val <- Q3 - Q1

    factor <- 5.0  # Very conservative
    lower_bound <- Q1 - factor * IQR_val
    upper_bound <- Q3 + factor * IQR_val

    outliers <- (values < lower_bound | values > upper_bound) & !is.na(values)

    return(list(
      outliers = outliers,
      method_used = "generic_conservative",
      statistics = list(
        total_outliers = sum(outliers, na.rm = TRUE),
        outlier_rate = mean(outliers, na.rm = TRUE)
      )
    ))
  }
}

#' Apply Basic Range Limits (Enhanced)
#'
#' Applies property-specific range constraints with enhanced property coverage.
#'
#' @param values Numeric vector
#' @param property_name Property name
#'
#' @return Constrained values
apply_basic_range_limits <- function(values, property_name) {

  if (all(is.na(values))) return(values)

  # Enhanced limits with comprehensive property coverage
  limits <- list(
    # Texture properties
    sandtotal = c(0, 100), claytotal = c(0, 100), silttotal = c(0, 100),
    sand_vc = c(0, 80), sand_c = c(0, 80), sand_m = c(0, 80),
    sand_f = c(0, 80), sand_vf = c(0, 80),

    # Physical properties
    dbovendry = c(0.3, 3.0), dbthirdbar = c(0.3, 3.0),
    partdensity = c(2.0, 3.2),

    # Water retention
    wthirdbar = c(0, 70), wfifteenbar = c(0, 60), awc = c(0, 0.5),

    # Chemical properties
    cec7 = c(0, 200), cec82 = c(0, 200), ecec = c(0, 200),
    ph1to1h2o = c(2.5, 11.0), ph01mcacl2 = c(2.0, 10.5),

    # Organic matter and carbon
    om = c(0, 100), oc = c(0, 60),

    # Rock fragments
    rfv = c(0, 95), fragvol = c(0, 95)
  )

  if (property_name %in% names(limits)) {
    range_limits <- limits[[property_name]]
    values[values < range_limits[1] | values > range_limits[2]] <- NA
  } else {
    # Generic non-negative constraint
    values[values < 0] <- NA
  }

  return(values)
}

#' Standardize Horizon Name
#'
#' Cleans and standardizes horizon names for better matching.
#'
#' @param hzname Raw horizon name
#'
#' @return Standardized horizon name
standardize_horizon_name <- function(hzname) {
  if (is.na(hzname)) return("")

  # Convert to string and clean
  hzname <- toupper(trimws(as.character(hzname)))

  # Remove numbers at the end
  hzname <- stringr::str_replace(hzname, "\\d+$", "")

  # Remove common punctuation but keep important ones
  hzname <- stringr::str_replace_all(hzname, "[^\\w/]", "")

  return(hzname)
}

#' Calculate Horizon Similarity
#'
#' Calculates similarity between two horizon names using pedological knowledge.
#'
#' @param hz1,hz2 Horizon names to compare
#'
#' @return Similarity score from 0 to 1
calculate_horizon_similarity <- function(hz1, hz2) {
  if (hz1 == hz2) return(1.0)
  if (hz1 == "" || hz2 == "") return(0.0)

  # Main horizon letter match
  main_hz1 <- ifelse(nchar(hz1) > 0, substr(hz1, 1, 1), "")
  main_hz2 <- ifelse(nchar(hz2) > 0, substr(hz2, 1, 1), "")

  if (main_hz1 == main_hz2) {
    base_score <- 0.8

    # Bonus for additional character matches
    chars1 <- strsplit(hz1, "")[[1]]
    chars2 <- strsplit(hz2, "")[[1]]
    common_chars <- intersect(chars1, chars2)
    bonus <- length(common_chars) / max(length(chars1), length(chars2)) * 0.2

    return(min(1.0, base_score + bonus))
  }

  # Related horizon groups
  related_groups <- list(
    c('A', 'AP', 'AE'),
    c('E', 'EB', 'BE'),
    c('B', 'BT', 'BW', 'BC', 'BS'),
    c('C', 'CB', 'CR'),
    c('O', 'OA', 'OE')
  )

  for (group in related_groups) {
    if (main_hz1 %in% group && main_hz2 %in% group) {
      return(0.6)
    }
  }

  return(0.0)
}

#' Calculate Property Lower Bound
#'
#' Calculates appropriate lower bound for property ranges.
#'
#' @param row Data frame row
#' @param property_name Property name
#' @param learned_ranges Learned ranges from data
#' @param context_ranges Contextual ranges
#' @param property_config Property configuration
#'
#' @return Lower bound value
calculate_property_lower_bound <- function(row, property_name, learned_ranges, context_ranges, property_config) {

  r_col <- paste0(property_name, "_r")
  r_value <- as.numeric(row[[r_col]])

  if (is.na(r_value) || !is.finite(r_value)) return(NA)

  # Priority-based spread estimation
  spread_estimates <- list()

  # 1. Horizon-specific learned range
  if ('hzname' %in% names(row) && !is.na(row[['hzname']])) {
    hz_key <- paste0("hzname_", row[['hzname']])
    if (hz_key %in% names(learned_ranges) && learned_ranges[[hz_key]]$sample_size >= 3) {
      spread <- learned_ranges[[hz_key]]$lower_spread_median
      if (!is.na(spread) && spread > 0) {
        spread_estimates <- append(spread_estimates, list(list(method = 'horizon_learned', spread = spread)))
      }
    }
  }

  # 2. Depth zone learned range
  if ('hzdepb_r' %in% names(row) && !is.na(row[['hzdepb_r']])) {
    depth <- row[['hzdepb_r']]
    zone_key <- if (depth <= 30) 'depth_surface' else if (depth <= 100) 'depth_subsurface' else 'depth_deep'

    if (zone_key %in% names(learned_ranges) && learned_ranges[[zone_key]]$sample_size >= 3) {
      spread <- learned_ranges[[zone_key]]$lower_spread_median
      if (!is.na(spread) && spread > 0) {
        spread_estimates <- append(spread_estimates, list(list(method = 'depth_learned', spread = spread)))
      }
    }
  }

  # 3. Overall dataset learned range
  if ('overall' %in% names(learned_ranges)) {
    spread <- learned_ranges[['overall']]$lower_spread_median
    if (!is.na(spread) && spread > 0) {
      spread_estimates <- append(spread_estimates, list(list(method = 'overall_learned', spread = spread)))
    }
  }

  # 4. Contextual knowledge
  ctx_spread <- get_contextual_spread(row, property_name, context_ranges, property_config, 'lower')
  if (!is.na(ctx_spread)) {
    spread_estimates <- append(spread_estimates, list(list(method = 'contextual', spread = ctx_spread)))
  }

  # 5. Fallback
  if (length(spread_estimates) == 0) {
    fallback <- property_config$fallback_range %||% 5
    spread_estimates <- append(spread_estimates, list(list(method = 'fallback', spread = fallback)))
  }

  # Use the first estimate
  spread <- spread_estimates[[1]]$spread
  if (!is.numeric(spread) || !is.finite(spread) || spread < 0) {
    spread <- property_config$fallback_range %||% 5
  }

  lower_bound <- r_value - spread

  # Apply constraints
  if (!is.null(property_config$typical_range)) {
    lower_bound <- max(lower_bound, property_config$typical_range[1])
  } else {
    lower_bound <- max(lower_bound, 0)
  }

  return(lower_bound)
}

#' Calculate Property Upper Bound
#'
#' Calculates appropriate upper bound for property ranges.
#'
#' @param row Data frame row
#' @param property_name Property name
#' @param learned_ranges Learned ranges from data
#' @param context_ranges Contextual ranges
#' @param property_config Property configuration
#'
#' @return Upper bound value
calculate_property_upper_bound <- function(row, property_name, learned_ranges, context_ranges, property_config) {

  r_col <- paste0(property_name, "_r")
  r_value <- as.numeric(row[[r_col]])

  if (is.na(r_value) || !is.finite(r_value)) return(NA)

  # Similar logic to lower bound but for upper spreads
  spread_estimates <- list()

  # Follow same priority order as lower bound
  if ('hzname' %in% names(row) && !is.na(row[['hzname']])) {
    hz_key <- paste0("hzname_", row[['hzname']])
    if (hz_key %in% names(learned_ranges) && learned_ranges[[hz_key]]$sample_size >= 3) {
      spread <- learned_ranges[[hz_key]]$upper_spread_median
      if (!is.na(spread) && spread > 0) {
        spread_estimates <- append(spread_estimates, list(list(method = 'horizon_learned', spread = spread)))
      }
    }
  }

  if ('hzdepb_r' %in% names(row) && !is.na(row[['hzdepb_r']])) {
    depth <- row[['hzdepb_r']]
    zone_key <- if (depth <= 30) 'depth_surface' else if (depth <= 100) 'depth_subsurface' else 'depth_deep'

    if (zone_key %in% names(learned_ranges) && learned_ranges[[zone_key]]$sample_size >= 3) {
      spread <- learned_ranges[[zone_key]]$upper_spread_median
      if (!is.na(spread) && spread > 0) {
        spread_estimates <- append(spread_estimates, list(list(method = 'depth_learned', spread = spread)))
      }
    }
  }

  if ('overall' %in% names(learned_ranges)) {
    spread <- learned_ranges[['overall']]$upper_spread_median
    if (!is.na(spread) && spread > 0) {
      spread_estimates <- append(spread_estimates, list(list(method = 'overall_learned', spread = spread)))
    }
  }

  ctx_spread <- get_contextual_spread(row, property_name, context_ranges, property_config, 'upper')
  if (!is.na(ctx_spread)) {
    spread_estimates <- append(spread_estimates, list(list(method = 'contextual', spread = ctx_spread)))
  }

  if (length(spread_estimates) == 0) {
    fallback <- property_config$fallback_range %||% 5
    spread_estimates <- append(spread_estimates, list(list(method = 'fallback', spread = fallback)))
  }

  spread <- spread_estimates[[1]]$spread
  if (!is.numeric(spread) || !is.finite(spread) || spread < 0) {
    spread <- property_config$fallback_range %||% 5
  }

  upper_bound <- r_value + spread

  # Apply constraints
  if (!is.null(property_config$typical_range)) {
    upper_bound <- min(upper_bound, property_config$typical_range[2])
  }

  return(upper_bound)
}

#' Get Contextual Spread
#'
#' Returns property-specific contextual spread based on soil characteristics.
#'
#' @param row Data frame row
#' @param property_name Property name
#' @param context_ranges Contextual ranges
#' @param property_config Property configuration
#' @param bound_type 'lower' or 'upper'
#'
#' @return Contextual spread value
get_contextual_spread <- function(row, property_name, context_ranges, property_config, bound_type) {

  # Property-specific spread calculation with pedological knowledge
  if (property_config$type == 'texture' && property_name %in% names(context_ranges)) {
    return(context_ranges[[property_name]]$typical_spread)
  }

  else if (property_config$type == 'bulk_density') {
    if ('hzname' %in% names(row) && !is.na(row[['hzname']])) {
      hzname <- toupper(as.character(row[['hzname']]))
      if (substr(hzname, 1, 1) == 'A' && 'A_horizons' %in% names(context_ranges)) {
        return(context_ranges[['A_horizons']]$typical_spread)
      } else if (substr(hzname, 1, 1) == 'B' && 'B_horizons' %in% names(context_ranges)) {
        return(context_ranges[['B_horizons']]$typical_spread)
      } else if (substr(hzname, 1, 1) == 'C' && 'C_horizons' %in% names(context_ranges)) {
        return(context_ranges[['C_horizons']]$typical_spread)
      }
    }
    return(property_config$fallback_range)
  }

  # Fallback to configuration default
  return(property_config$fallback_range %||% 5)
}

#' Impute Rock Fragment Volume (RFV) Values for One Row
#'
#' Row-level RFV imputation: texture-informed defaults when `rfv_r` is
#' missing, floor handling for near-zero values, and a simple
#' proportional-spread `_l/_h` estimate (+/-30%, clamped to `[0.05, 85]`)
#' for valid representative values.
#'
#' @param row A one-row data frame or list with `claytotal_r`/`sandtotal_r`/
#'   `silttotal_r`/`rfv_r` fields.
#' @return A one-row data frame with `rfv_l`/`rfv_r`/`rfv_h` set.
#' @export
impute_rfv_values <- function(row) {
  if (is.data.frame(row)) row <- as.list(row[1, ])

  # Get current values
  claytotal_r <- as.numeric(row[["claytotal_r"]])
  sandtotal_r <- as.numeric(row[["sandtotal_r"]])
  silttotal_r <- as.numeric(row[["silttotal_r"]])
  rfv_r <- as.numeric(row[["rfv_r"]])

  # Handle missing RFV
  if (is.na(rfv_r)) {
    # Use texture to inform estimate
    if (!is.na(claytotal_r) && !is.na(silttotal_r) && !is.na(sandtotal_r)) {
      texture_sum <- sandtotal_r + silttotal_r + claytotal_r
      if (texture_sum >= 95 && texture_sum <= 105) {
        row[["rfv_r"]] <- 0.5
        row[["rfv_l"]] <- 0.1
        row[["rfv_h"]] <- 1.0
        return(as.data.frame(row))
      }
    }

    # Default estimate
    row[["rfv_r"]] <- 2.0
    row[["rfv_l"]] <- 0.5
    row[["rfv_h"]] <- 5.0
    return(as.data.frame(row))
  }

  # Handle zero/low values
  if (rfv_r <= 0.01) {
    row[["rfv_r"]] <- 0.1
    row[["rfv_l"]] <- 0.05
    row[["rfv_h"]] <- 0.5
    return(as.data.frame(row))
  }

  # Calculate ranges for valid values
  rfv_r <- max(0.1, min(rfv_r, 85))
  row[["rfv_r"]] <- rfv_r

  # Simple range calculation
  spread <- rfv_r * 0.3
  row[["rfv_l"]] <- max(0.05, rfv_r - spread)
  row[["rfv_h"]] <- min(85, rfv_r + spread)

  return(as.data.frame(row))
}

#' Saxton-Rawls Water Retention Pedotransfer Function (Single Horizon)
#'
#' Estimates field capacity and wilting point water content from texture,
#' bulk density, rock fragment volume, and organic matter, via the
#' Saxton-Rawls pedotransfer equations (simplified/RFV-corrected variant).
#' Real, checkable math (not a placeholder) - inputs are clamped to
#' physically plausible ranges and texture percentages are renormalized to
#' sum to 100 when off by more than 5 points.
#'
#' @param sand_pct,clay_pct,silt_pct Texture percentages (0-100).
#' @param bulk_density Bulk density (g/cm^3), clamped to `[0.6, 2.5]`.
#' @param rfv_pct Rock fragment volume percentage, clamped to `[0, 95]`.
#' @param om_pct Organic matter percentage, clamped to `[0.1, 50]`.
#' @return A list of estimated water retention values (field capacity,
#'   wilting point, and their `_l/_h` spread) - see the function body for
#'   the exact returned fields.
#' @export
calculate_saxton_rawls_single <- function(sand_pct, clay_pct, silt_pct, bulk_density, rfv_pct = 0, om_pct = 2) {

  # Input validation
  sand_pct <- max(0, min(sand_pct, 100))
  clay_pct <- max(0, min(clay_pct, 100))
  silt_pct <- max(0, min(silt_pct, 100))
  bulk_density <- max(0.6, min(bulk_density, 2.5))
  rfv_pct <- max(0, min(rfv_pct, 95))
  om_pct <- max(0.1, min(om_pct, 50))

  # Normalize texture
  texture_sum <- sand_pct + clay_pct + silt_pct
  if (abs(texture_sum - 100) > 5) {
    sand_pct <- sand_pct / texture_sum * 100
    clay_pct <- clay_pct / texture_sum * 100
    silt_pct <- silt_pct / texture_sum * 100
  }

  # Convert to fractions
  sand <- sand_pct / 100
  clay <- clay_pct / 100
  om <- om_pct / 100

  # Saxton-Rawls equations (simplified)
  theta_s <- -0.251 * sand + 0.195 * clay + 0.011 * om + 0.006 * sand * om -
    0.027 * clay * om + 0.452 * sand * clay + 0.299

  # Field capacity
  A <- exp(log(33) + 1.54 * sand + 0.95 * clay + 0.025 * om - 0.351 * sand * om -
             0.023 * clay * om - 0.427 * sand * clay + 0.015 * sand^2 * clay^2)
  fc_gravimetric <- theta_s * (A / (A + 0.31))^0.17

  # Wilting point
  B <- exp(log(1500) + 0.02 * clay^2 + 0.14 * sand - 0.0002 * sand^2 * clay -
             0.002 * clay^2 * sand - 0.0002 * clay^2 * om + 0.0003 * clay^2 * sand * om)
  wp_gravimetric <- theta_s * (B / (B + 0.31))^0.17

  # Convert to volumetric and apply RFV correction
  rfv_fraction <- rfv_pct / 100
  fc_volumetric <- fc_gravimetric * bulk_density * (1 - rfv_fraction) * 100
  wp_volumetric <- wp_gravimetric * bulk_density * (1 - rfv_fraction) * 100

  # Apply constraints
  field_capacity <- max(3, min(fc_volumetric, 65))
  wilting_point <- max(1, min(wp_volumetric, 45))

  if (wilting_point >= field_capacity) {
    wilting_point <- field_capacity * 0.6
  }

  # Calculate ranges (+/-15% uncertainty)
  fc_spread <- field_capacity * 0.15
  wp_spread <- wilting_point * 0.15

  return(list(
    field_capacity = round(field_capacity, 2),
    wilting_point = round(wilting_point, 2),
    field_capacity_l = round(max(2, field_capacity - fc_spread), 2),
    field_capacity_h = round(min(65, field_capacity + fc_spread), 2),
    wilting_point_l = round(max(1, wilting_point - wp_spread), 2),
    wilting_point_h = round(min(45, wilting_point + wp_spread), 2),
    available_water_capacity = round(field_capacity - wilting_point, 2)
  ))
}

#' Depth Weighted Property Infill
#'
#' @param group Data frame group
#' @param property_col Property column name
#' @param problematic_mask Logical mask for problematic horizons
#'
#' @return Data frame with infilled values
depth_weighted_property_infill <- function(group, property_col, problematic_mask) {

  if (!all(c("hzdept_r", "hzdepb_r") %in% names(group))) return(group)

  # Validate mask length
  if (length(problematic_mask) != nrow(group)) {
    warning("Mask length mismatch in depth_weighted_property_infill")
    return(group)
  }

  suitable_mask <- if ("unsuitable_horizon" %in% names(group)) {
    !group$unsuitable_horizon
  } else {
    rep(TRUE, nrow(group))
  }

  problematic_indices <- which(problematic_mask & is.na(group[[property_col]]))
  source_indices <- which(suitable_mask & !is.na(group[[property_col]]))

  if (length(source_indices) == 0) return(group)

  for (idx in problematic_indices) {
    # Calculate target midpoint
    target_mid <- (group$hzdept_r[idx] + group$hzdepb_r[idx]) / 2

    # Skip if target depth is invalid
    if (is.na(target_mid)) next

    # Find sources within 20cm depth tolerance
    source_mids <- (group$hzdept_r[source_indices] + group$hzdepb_r[source_indices]) / 2
    depth_diffs <- abs(target_mid - source_mids)

    # Handle NA values in depth differences
    valid_diffs <- !is.na(depth_diffs)
    if (!any(valid_diffs)) next

    depth_diffs <- depth_diffs[valid_diffs]
    valid_source_indices <- source_indices[valid_diffs]

    within_tolerance <- depth_diffs <= 20

    # FIX: Check if any are within tolerance and handle NA values properly
    if (length(within_tolerance) > 0 && any(within_tolerance, na.rm = TRUE)) {
      matching_indices <- valid_source_indices[within_tolerance]
      matching_values <- group[[property_col]][matching_indices]
      matching_diffs <- depth_diffs[within_tolerance]

      # Remove any NA values from matching data
      valid_matches <- !is.na(matching_values) & !is.na(matching_diffs)
      if (any(valid_matches)) {
        final_values <- matching_values[valid_matches]
        final_diffs <- matching_diffs[valid_matches]

        weights <- 1 / (1 + final_diffs)
        group[[property_col]][idx] <- weighted.mean(final_values, weights)

        # Update tracking
        if (!"infill_method" %in% names(group)) {
          group$infill_method <- ""
        }

        method_info <- paste0(property_col, ":depth_weighted(n=", length(final_values), "); ")
        if (is.na(group$infill_method[idx]) || group$infill_method[idx] == "") {
          group$infill_method[idx] <- method_info
        } else {
          group$infill_method[idx] <- paste0(group$infill_method[idx], method_info)
        }
      }
    }
  }

  return(group)
}

#' Within Component Interpolation
#'
#' @param group Data frame group
#' @param property_col Property column name
#'
#' @return Data frame with interpolated values
within_component_property_interpolation <- function(group, property_col) {

  if (!all(c("hzdept_r", "hzdepb_r") %in% names(group))) return(group)

  suitable_mask <- if ("unsuitable_horizon" %in% names(group)) {
    !group$unsuitable_horizon
  } else {
    rep(TRUE, nrow(group))
  }

  suitable_data <- group[suitable_mask, ]

  if (nrow(suitable_data) < 2) return(group)

  valid_data_mask <- !is.na(suitable_data[[property_col]]) &
    !is.na(suitable_data$hzdepb_r)
  missing_data_mask <- is.na(suitable_data[[property_col]]) &
    !is.na(suitable_data$hzdepb_r)

  if (sum(valid_data_mask) >= 2 && any(missing_data_mask)) {
    valid_indices <- which(valid_data_mask)
    missing_indices <- which(missing_data_mask)

    tryCatch({
      interpolated_values <- approx(
        x = suitable_data$hzdepb_r[valid_indices],
        y = suitable_data[[property_col]][valid_indices],
        xout = suitable_data$hzdepb_r[missing_indices],
        rule = 2
      )$y

      # Update only the suitable data
      suitable_data[[property_col]][missing_indices] <- interpolated_values

      # Update tracking for interpolated values
      if (!"infill_method" %in% names(suitable_data)) {
        suitable_data$infill_method <- ""
      }

      method_info <- paste0(property_col, ":interpolation; ")
      suitable_data$infill_method[missing_indices] <- paste0(
        suitable_data$infill_method[missing_indices], method_info)

      # Put suitable data back into the full group
      group[suitable_mask, ] <- suitable_data

    }, error = function(e) {
      warning(paste("Interpolation failed for", property_col, ":", e$message))
    })
  }

  return(group)
}

#' Cross-Component Property Interpolation
#'
#' Fills missing values in suitable horizons using data from other soil
#' components (other `cokey` groups within the same `group`) at similar
#' depths, weighted by depth proximity. Only suitable horizons are used as
#' sources.
#'
#' @param group Data frame group to process
#' @param property_col Name of the property column to infill
#'
#' @return Data frame with infilled values
cross_component_property_interpolation <- function(group, property_col) {

  if (!is.data.frame(group) || nrow(group) == 0) return(group)
  if (!property_col %in% names(group)) return(group)

  required_depth_cols <- c("hzdept_r", "hzdepb_r")
  if (!all(required_depth_cols %in% names(group))) return(group)

  if (!"infill_method" %in% names(group)) {
    group$infill_method <- ""
  }

  suitable_mask <- if ("unsuitable_horizon" %in% names(group)) {
    !group$unsuitable_horizon
  } else {
    rep(TRUE, nrow(group))
  }

  missing_mask <- is.na(group[[property_col]]) &
    suitable_mask &
    !is.na(group$hzdept_r) & !is.na(group$hzdepb_r) &
    is.finite(group$hzdept_r) & is.finite(group$hzdepb_r)

  if (!any(missing_mask)) return(group)

  source_mask <- suitable_mask &
    !is.na(group[[property_col]]) &
    !is.na(group$hzdept_r) & !is.na(group$hzdepb_r) &
    is.finite(group[[property_col]]) &
    is.finite(group$hzdept_r) & is.finite(group$hzdepb_r)

  source_indices <- which(source_mask)
  if (length(source_indices) == 0) return(group)

  source_tops <- group$hzdept_r[source_indices]
  source_bottoms <- group$hzdepb_r[source_indices]
  source_mids <- (source_tops + source_bottoms) / 2
  source_values <- group[[property_col]][source_indices]

  valid_source_mask <- is.finite(source_mids)
  if (!any(valid_source_mask)) return(group)

  source_indices <- source_indices[valid_source_mask]
  source_mids <- source_mids[valid_source_mask]
  source_values <- source_values[valid_source_mask]

  missing_indices <- which(missing_mask)
  depth_tolerance <- 15  # 15cm tolerance for similar depths, matching legacy reference behavior

  for (idx in missing_indices) {
    target_top <- group$hzdept_r[idx]
    target_bottom <- group$hzdepb_r[idx]

    if (!is.finite(target_top) || !is.finite(target_bottom) || target_top >= target_bottom) next

    target_mid <- (target_top + target_bottom) / 2
    if (!is.finite(target_mid)) next

    depth_diffs <- abs(target_mid - source_mids)
    within_tolerance <- (depth_diffs <= depth_tolerance) & (source_indices != idx)

    if (!any(within_tolerance)) next

    matching_values <- source_values[within_tolerance]
    matching_diffs <- depth_diffs[within_tolerance]
    weights <- 1 / (1 + matching_diffs)

    valid_match_mask <- is.finite(weights) & is.finite(matching_values) & weights > 0
    if (!any(valid_match_mask)) next

    final_values <- matching_values[valid_match_mask]
    final_weights <- weights[valid_match_mask]

    weighted_value <- if (length(final_values) == 1) {
      final_values[1]
    } else {
      stats::weighted.mean(final_values, final_weights)
    }

    if (is.finite(weighted_value)) {
      group[[property_col]][idx] <- weighted_value

      depth_info <- paste0(round(target_mid, 1), "cm")
      method_info <- paste0(property_col, ":cross_comp(", depth_info, ",n=",
                            length(final_values), ",tol=", depth_tolerance, "cm); ")

      group$infill_method[idx] <- paste0(group$infill_method[idx], method_info)
    }
  }

  return(group)
}

#' Related Property Estimation
#'
#' Estimates missing values from pedologically-related properties, using
#' property-type-specific relationships (e.g. texture sum-to-100 constraint,
#' clay/organic-matter based CEC estimation, clay-based water retention,
#' horizon/organic-matter adjusted pH, depth/horizon/clay adjusted organic
#' matter, texture-adjusted bulk density). Only suitable horizons are
#' estimated.
#'
#' @param group Data frame group to process
#' @param property_name Name of the property to estimate
#' @param property_config Property configuration list
#'
#' @return Data frame with estimated values
related_property_estimation <- function(group, property_name, property_config) {

  property_col <- paste0(property_name, "_r")
  if (!property_col %in% names(group)) return(group)

  suitable_mask <- if ("unsuitable_horizon" %in% names(group)) {
    !group$unsuitable_horizon
  } else {
    rep(TRUE, nrow(group))
  }

  if (is.null(property_config$related_properties)) return(group)

  missing_mask <- is.na(group[[property_col]]) & suitable_mask
  if (!any(missing_mask)) return(group)

  if (!"infill_method" %in% names(group)) {
    group$infill_method <- ""
  }

  mark_estimated <- function(group, idx, tag) {
    group$infill_method[idx] <- paste0(group$infill_method[idx], property_col, ":", tag, "; ")
    group
  }

  # Texture properties - use sum constraint (sand + clay + silt ~= 100%)
  if (property_config$type == 'texture') {
    texture_cols <- paste0(property_config$related_properties, "_r")
    available_cols <- intersect(texture_cols, names(group))

    if (length(available_cols) >= 2) {
      for (idx in which(missing_mask)) {
        other_values <- unlist(group[idx, available_cols])
        if (sum(!is.na(other_values)) >= 2) {
          sum_others <- sum(other_values, na.rm = TRUE)
          estimated_value <- max(0, min(100, 100 - sum_others))
          group[[property_col]][idx] <- estimated_value
          group <- mark_estimated(group, idx, "related_texture_sum")
        }
      }
    }
  }

  # Water retention - use clay relationship if available
  else if (property_config$type == 'water_retention' && 'claytotal_r' %in% names(group)) {
    for (idx in which(missing_mask)) {
      clay_content <- group[['claytotal_r']][idx]
      if (!is.na(clay_content)) {
        estimated_value <- if (property_name == 'wthirdbar') {
          max(0, min(60, 0.3 * clay_content + 10))
        } else if (property_name == 'wfifteenbar') {
          max(0, min(40, 0.4 * clay_content + 2))
        } else {
          NA_real_
        }
        if (!is.na(estimated_value)) {
          group[[property_col]][idx] <- estimated_value
          group <- mark_estimated(group, idx, "related_clay")
        }
      }
    }
  }

  # CEC - use clay and organic matter relationships
  else if (property_config$type == 'cec') {
    for (idx in which(missing_mask)) {
      clay <- if ('claytotal_r' %in% names(group)) group$claytotal_r[idx] else NA_real_
      om <- if ('om_r' %in% names(group)) group$om_r[idx] else NA_real_

      if (!is.na(clay) || !is.na(om)) {
        estimated_cec <- 0
        if (!is.na(clay)) estimated_cec <- estimated_cec + (clay * 0.5)
        if (!is.na(om)) estimated_cec <- estimated_cec + (om * 20)
        estimated_cec <- max(2, estimated_cec)

        group[[property_col]][idx] <- max(0, min(100, estimated_cec))
        group <- mark_estimated(group, idx, "related_clay_om")
      }
    }
  }

  # pH - use horizon and organic matter context
  else if (property_config$type == 'ph') {
    for (idx in which(missing_mask)) {
      estimated_ph <- 6.2

      if ('hzname' %in% names(group) && !is.na(group$hzname[idx])) {
        hzname <- toupper(as.character(group$hzname[idx]))
        if (substr(hzname, 1, 1) == 'A') {
          estimated_ph <- estimated_ph - 0.3
        } else if (substr(hzname, 1, 1) == 'C') {
          estimated_ph <- estimated_ph + 0.2
        }
      }

      if ('om_r' %in% names(group) && !is.na(group$om_r[idx])) {
        if (group$om_r[idx] > 5) {
          estimated_ph <- estimated_ph - 0.4
        }
      }

      group[[property_col]][idx] <- max(3.0, min(10.0, estimated_ph))
      group <- mark_estimated(group, idx, "related_horizon_om")
    }
  }

  # Organic matter - use depth and horizon relationships
  else if (property_config$type == 'organic_matter') {
    for (idx in which(missing_mask)) {
      estimated_om <- 2.0

      if ('hzdept_r' %in% names(group) && !is.na(group$hzdept_r[idx])) {
        depth <- group$hzdept_r[idx]
        estimated_om <- if (depth <= 15) {
          3.5
        } else if (depth <= 30) {
          1.8
        } else if (depth <= 50) {
          0.8
        } else {
          0.3
        }
      }

      if ('hzname' %in% names(group) && !is.na(group$hzname[idx])) {
        hzname <- toupper(as.character(group$hzname[idx]))
        if (substr(hzname, 1, 1) == 'A') {
          estimated_om <- estimated_om * 1.5
        } else if (substr(hzname, 1, 1) == 'C') {
          estimated_om <- 0.2
        }
      }

      if ('claytotal_r' %in% names(group) && !is.na(group$claytotal_r[idx])) {
        clay_content <- group$claytotal_r[idx]
        if (clay_content > 35) {
          estimated_om <- estimated_om * 1.3
        } else if (clay_content < 15) {
          estimated_om <- estimated_om * 0.7
        }
      }

      group[[property_col]][idx] <- max(0, min(50, estimated_om))
      group <- mark_estimated(group, idx, "related_depth_horizon")
    }
  }

  # Bulk density - use texture relationship if available
  else if (property_config$type == 'bulk_density' &&
           any(c('sandtotal_r', 'claytotal_r') %in% names(group))) {
    for (idx in which(missing_mask)) {
      clay <- if ('claytotal_r' %in% names(group)) group$claytotal_r[idx] else NA_real_
      sand <- if ('sandtotal_r' %in% names(group)) group$sandtotal_r[idx] else NA_real_

      if (!is.na(clay) || !is.na(sand)) {
        base_bd <- 1.4
        if (!is.na(clay)) base_bd <- base_bd - (clay - 20) * 0.01
        if (!is.na(sand)) base_bd <- base_bd + (sand - 50) * 0.005

        group[[property_col]][idx] <- max(0.8, min(2.2, base_bd))
        group <- mark_estimated(group, idx, "related_texture")
      }
    }
  }

  return(group)
}

#' Calculate Depth Weighted Mean
#'
#' Handles NA values properly
#'
#' @param group Data frame group
#' @param property_col Property column name
#'
#' @return Depth-weighted mean value
calculate_depth_weighted_mean <- function(group, property_col) {

  suitable_mask <- if ("unsuitable_horizon" %in% names(group)) {
    !group$unsuitable_horizon
  } else {
    rep(TRUE, nrow(group))
  }

  valid_data <- group[!is.na(group[[property_col]]) & suitable_mask, ]

  if (nrow(valid_data) == 0) return(NA)

  if (all(c('hzdept_r', 'hzdepb_r') %in% names(valid_data))) {
    thicknesses <- valid_data$hzdepb_r - valid_data$hzdept_r

    # Remove invalid thicknesses
    valid_thickness_mask <- !is.na(thicknesses) & thicknesses > 0

    if (any(valid_thickness_mask)) {
      valid_thicknesses <- thicknesses[valid_thickness_mask]
      valid_values <- valid_data[[property_col]][valid_thickness_mask]

      if (sum(valid_thicknesses) > 0) {
        return(weighted.mean(valid_values, valid_thicknesses))
      }
    }
  }

  return(mean(valid_data[[property_col]], na.rm = TRUE))
}

#' Process Property Group Function
#'
#' @param group Data frame group to process
#' @param property_name Name of the property to infill
#' @param problematic_mask Logical vector for problematic horizons (IGNORED - recalculated)
#' @param property_config Property configuration
#' @param max_depth Maximum depth for processing
#' @param verbose Whether to provide progress messages
#'
#' @return Data frame with infilled property data
process_property_group <- function(group, property_name, problematic_mask,
                                   property_config, max_depth, verbose) {

  # RECALCULATE problematic_mask for this group to avoid length mismatch
  property_col <- paste0(property_name, "_r")

  # Ensure unsuitable horizon detection for this group
  if (!"unsuitable_horizon" %in% names(group)) {
    group$unsuitable_horizon <- is_unsuitable(group, hzname_col = "hzname")
  }

  # Apply depth constraints for this group
  if ("hzdepb_r" %in% names(group)) {
    depth_mask <- group$hzdepb_r <= max_depth
  } else {
    depth_mask <- rep(TRUE, nrow(group))
  }

  # RECALCULATE problematic mask for this specific group
  group_problematic_mask <- is.na(group[[property_col]]) &
    !group$unsuitable_horizon &
    depth_mask

  # Apply recovery strategies with corrected mask
  group <- infill_missing_property_data(group, property_name,
                                        group_problematic_mask, property_config)

  # Set completion flag
  still_missing <- is.na(group[[property_col]])
  group$property_data_complete <- !(still_missing & !group$unsuitable_horizon)

  return(group)
}

#' Process Property Group Fallback Function
#'
#' Companion to [process_property_group()]: applies `apply_group_fallback_mean()`
#' (Strategy 6, last resort) per grouping unit, after the whole-dataset
#' Strategy 4/5 passes (cross-component interpolation, related-property
#' estimation) have already run in [infill_soil_property()]. Recalculates its
#' own problematic mask the same way [process_property_group()] does, so it
#' only fills cells that are still genuinely missing at this point.
#'
#' @param group Data frame group to process
#' @param property_name Name of the property to infill
#' @param property_config Property configuration
#' @param max_depth Maximum depth for processing
#'
#' @return Data frame with fallback-filled property data
process_property_group_fallback <- function(group, property_name, property_config, max_depth) {

  property_col <- paste0(property_name, "_r")

  if (!"unsuitable_horizon" %in% names(group)) {
    group$unsuitable_horizon <- is_unsuitable(group, hzname_col = "hzname")
  }

  if ("hzdepb_r" %in% names(group)) {
    depth_mask <- group$hzdepb_r <= max_depth
  } else {
    depth_mask <- rep(TRUE, nrow(group))
  }

  group_problematic_mask <- is.na(group[[property_col]]) &
    !group$unsuitable_horizon &
    depth_mask

  group <- apply_group_fallback_mean(group, property_col, group_problematic_mask)

  still_missing <- is.na(group[[property_col]])
  group$property_data_complete <- !(still_missing & !group$unsuitable_horizon)

  return(group)
}
