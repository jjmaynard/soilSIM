# System Utilities & Helpers - Module 8 (INTEGRATED VERSION)
# Core utility functions supporting all soil simulation workflow modules
# Provides data validation, transformations, I/O, logging, configuration management
# This module serves as the central utility hub for Modules 1 and 2

# ============================================================================
# 1. CORE DATA VALIDATION & QUALITY CONTROL FUNCTIONS
# ============================================================================

#' Identify Unsuitable Horizons
#'
#' Identifies horizons that should be excluded from soil property modeling
#' based on horizon designation and characteristics. This is the centralized
#' function used by all modules.
#'
#' @param data Soil horizon data with horizon designation information
#' @param hzname_col Column name containing horizon designations (default = "hzname")
#' @param strict_mode Whether to use strict unsuitable criteria (default = TRUE)
#' @param custom_exclusions Additional horizon patterns to exclude
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's log level so \code{INFO}-level progress messages print for the duration of this call (default \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @return Logical vector indicating unsuitable horizons
#' @export
is_unsuitable <- function(data,
                          hzname_col = "hzname",
                          strict_mode = TRUE,
                          custom_exclusions = NULL,
                          verbose = getOption("ssurgo.verbose", FALSE)) {
  .old_log_cfg <- set_verbose_logging(verbose)
  on.exit(options(soil_workflow_log_config = .old_log_cfg), add = TRUE)

  if (!hzname_col %in% names(data)) {
    log_message("WARN", paste("Horizon name column '", hzname_col, "' not found. Returning all FALSE."), category = "Validation")
    return(rep(FALSE, nrow(data)))
  }

  hznames <- data[[hzname_col]]

  # Convert to character and handle NAs
  hznames <- as.character(hznames)
  hznames[is.na(hznames)] <- ""

  # Check what additional columns are available for enhanced detection
  has_desgnmaster <- "desgnmaster" %in% names(data)
  has_restriction_top <- "restriction_top" %in% names(data)
  has_hzdept_r <- "hzdept_r" %in% names(data)
  has_is_potentially_restrictive <- "is_potentially_restrictive" %in% names(data)
  has_missing_required_properties <- "missing_required_properties" %in% names(data)
  has_organic_horizon <- "organic_horizon" %in% names(data)

  n <- length(hznames)
  unsuitable <- rep(FALSE, n)

  # PRIORITY 1: Always check horizon names first (R, Cr, O should ALWAYS be unsuitable)
  for (i in seq_len(n)) {
    hz <- toupper(trimws(hznames[i]))
    master <- if (has_desgnmaster) toupper(trimws(as.character(data$desgnmaster[i]))) else NA

    if (is.na(hz) || hz == "") {
      if (!is.na(master) && master %in% c("R", "O")) {
        unsuitable[i] <- TRUE
      }
      next
    }

    # Match O, R, Cr, Cd, Cx using hzname or master designation
    if (!is.na(master) && master %in% c("R", "O")) {
      unsuitable[i] <- TRUE
      next
    }

    # Match R, Cr, Cd, Cx patterns (with optional numeric prefix)
    if (grepl("(^|[0-9]+)(R$|Cr|Cd|Cx)", hz)) {
      unsuitable[i] <- TRUE
      next
    }
    # Match restrictive descriptors
    if (grepl("\\bfragipan\\b|\\bduripan\\b|\\bpetrocalcic\\b|\\bpetrogypsic\\b|\\bcemented\\b|\\bindurated\\b", hz, ignore.case = TRUE)) {
      unsuitable[i] <- TRUE
      next
    }
    # Match cemented suffix 'm'
    if (nchar(hz) > 1 && grepl("m$", hz)) {
      unsuitable[i] <- TRUE
      next
    }
    # Match organic horizons
    if (grepl("(^|[0-9]+)O", hz)) {
      unsuitable[i] <- TRUE
      next
    }
  }

  # PRIORITY 2: Add high-detail restrictions if available
  if (has_is_potentially_restrictive && has_missing_required_properties && has_hzdept_r && has_restriction_top) {
    hz_below_restriction <- !is.na(data$hzdept_r) & !is.na(data$restriction_top) & data$hzdept_r >= data$restriction_top
    restriction_logic <- data$is_potentially_restrictive & data$missing_required_properties & hz_below_restriction
    unsuitable <- unsuitable | restriction_logic
  }

  # PRIORITY 3: Add organic horizon flags if available
  if (has_organic_horizon) {
    organic_logic <- data$organic_horizon
    unsuitable <- unsuitable | organic_logic
  }

  # Additional checks for specific conditions
  if (strict_mode) {
    # Check for very deep horizons that might be unsuitable
    if ("hzdept_r" %in% names(data)) {
      very_deep <- data$hzdept_r > 250  # Deeper than 250 cm
      unsuitable <- unsuitable | (very_deep & !is.na(very_deep))
    }

    # Additional patterns for strict mode
    strict_patterns <- c("^C$", "^C[0-9]", "BC", "CB", "2C", "3C", "4C")
    for (pattern in strict_patterns) {
      pattern_matches <- grepl(pattern, hznames, ignore.case = FALSE)
      unsuitable <- unsuitable | pattern_matches
    }
  }

  # Add custom exclusions
  if (!is.null(custom_exclusions)) {
    for (pattern in custom_exclusions) {
      custom_matches <- grepl(pattern, hznames, ignore.case = FALSE)
      unsuitable <- unsuitable | custom_matches
    }
  }

  # Log unsuitable horizon statistics
  n_unsuitable <- sum(unsuitable, na.rm = TRUE)
  pct_unsuitable <- round(n_unsuitable / length(unsuitable) * 100, 1)

  if (n_unsuitable > 0) {
    log_message("INFO", paste("Identified", n_unsuitable, "unsuitable horizons (", pct_unsuitable, "%)"), category = "Validation")

    # Show unsuitable patterns found
    unsuitable_hznames <- unique(hznames[unsuitable])
    if (length(unsuitable_hznames) <= 10) {
      log_message("DEBUG", paste("Unsuitable horizon types:", paste(unsuitable_hznames, collapse = ", ")), category = "Validation")
    } else {
      log_message("DEBUG", paste("Unsuitable horizon types:", paste(head(unsuitable_hznames, 10), collapse = ", "), "..."), category = "Validation")
    }
  }

  return(unsuitable)
}

#' Validate Data Quality
#'
#' Comprehensive data quality assessment for soil data including missing values,
#' outliers, and data consistency checks.
#'
#' @param data Input soil data
#' @param required_columns Required columns for the analysis
#' @param numeric_columns Columns that should be numeric
#' @param quality_thresholds Quality assessment thresholds
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's log level so \code{INFO}-level progress messages print for the duration of this call (default \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @return Data quality assessment results
#' @export
validate_data_quality <- function(data,
                                  required_columns = character(0),
                                  numeric_columns = character(0),
                                  quality_thresholds = NULL,
                                  verbose = getOption("ssurgo.verbose", FALSE)) {
  .old_log_cfg <- set_verbose_logging(verbose)
  on.exit(options(soil_workflow_log_config = .old_log_cfg), add = TRUE)

  log_message("INFO", "=== DATA QUALITY VALIDATION ===", category = "Validation")
  log_message("INFO", paste("Rows:", nrow(data), ", Columns:", ncol(data)), category = "Validation")

  # Initialize validation structure
  validation_result <- list(
    overall_quality = list(),
    column_assessments = list(),
    missing_data_summary = list(),
    numeric_data_summary = list(),
    outlier_summary = list(),
    consistency_summary = list()
  )

  # FIXED: Safe quality thresholds with comprehensive defaults
  if (is.null(quality_thresholds)) {
    quality_thresholds <- list(
      max_missing_percentage = 0.2,  # 20% max missing
      min_completeness = 0.8,        # 80% min completeness
      max_outlier_percentage = 0.05, # 5% max outliers
      min_numeric_columns = 1        # At least 1 numeric column
    )
  } else {
    # Ensure all required thresholds exist
    quality_thresholds$max_missing_percentage <- quality_thresholds$max_missing_percentage %||% 0.2
    quality_thresholds$min_completeness <- quality_thresholds$min_completeness %||% 0.8
    quality_thresholds$max_outlier_percentage <- quality_thresholds$max_outlier_percentage %||% 0.05
    quality_thresholds$min_numeric_columns <- quality_thresholds$min_numeric_columns %||% 1
  }

  # Step 1: Missing data assessment with FIXED mutate logic
  log_message("DEBUG", "Missing data assessment", category = "Validation")

  tryCatch({
    # FIXED: Safe missing data assessment
    column_stats <- data.frame(
      column_name = names(data),
      stringsAsFactors = FALSE
    )

    # Calculate missing percentages safely
    missing_percentages <- sapply(data, function(col) {
      if (length(col) == 0) return(1.0)  # Empty column = 100% missing
      mean(is.na(col))
    })

    # Ensure we have the right length
    if (length(missing_percentages) != nrow(column_stats)) {
      log_message("WARN", "Missing percentage calculation mismatch, using safe defaults", category = "Validation")
      missing_percentages <- rep(0, nrow(column_stats))
    }

    column_stats$missing_percentage <- missing_percentages

    # FIXED: Safe quality flag creation with proper vector operations
    threshold_value <- quality_thresholds$max_missing_percentage

    # Ensure threshold is a single value
    if (length(threshold_value) == 0) {
      threshold_value <- 0.2  # Default
    } else if (length(threshold_value) > 1) {
      threshold_value <- threshold_value[1]  # Take first value
    }

    column_stats$quality_flag <- column_stats$missing_percentage > threshold_value

    # Count problematic columns
    high_missing_cols <- sum(column_stats$quality_flag, na.rm = TRUE)

    log_message("DEBUG", paste("Columns with >", round(threshold_value * 100), "% missing:", high_missing_cols), category = "DataAssessment")

    validation_result$missing_data_summary <- list(
      total_columns = ncol(data),
      high_missing_columns = high_missing_cols,
      avg_missing_percentage = mean(column_stats$missing_percentage, na.rm = TRUE),
      max_missing_percentage = max(column_stats$missing_percentage, na.rm = TRUE),
      column_details = column_stats
    )

  }, error = function(e) {
    log_message("WARN", paste("Missing data assessment error:", e$message), category = "Validation")

    # Fallback missing data summary
    validation_result$missing_data_summary <<- list(
      total_columns = ncol(data),
      high_missing_columns = 0,
      avg_missing_percentage = 0,
      max_missing_percentage = 0,
      assessment_failed = TRUE,
      error_message = e$message
    )
  })

  # Step 2: Numeric data validation with enhanced error handling
  log_message("DEBUG", "Numeric data validation", category = "Validation")

  tryCatch({
    if (length(numeric_columns) > 0) {
      # Filter to existing numeric columns
      existing_numeric <- intersect(numeric_columns, names(data))

      if (length(existing_numeric) > 0) {
        numeric_stats <- list()

        for (col in existing_numeric) {
          col_data <- data[[col]]
          finite_data <- col_data[is.finite(col_data)]

          numeric_stats[[col]] <- list(
            total_values = length(col_data),
            finite_values = length(finite_data),
            missing_values = sum(is.na(col_data)),
            infinite_values = sum(!is.finite(col_data) & !is.na(col_data)),
            completeness = length(finite_data) / length(col_data)
          )
        }

        validation_result$numeric_data_summary <- list(
          total_numeric_columns = length(existing_numeric),
          requested_numeric_columns = length(numeric_columns),
          column_stats = numeric_stats,
          avg_completeness = mean(sapply(numeric_stats, function(x) x$completeness), na.rm = TRUE)
        )
      } else {
        validation_result$numeric_data_summary <- list(
          total_numeric_columns = 0,
          requested_numeric_columns = length(numeric_columns),
          no_matching_columns = TRUE
        )
      }
    } else {
      validation_result$numeric_data_summary <- list(
        total_numeric_columns = 0,
        requested_numeric_columns = 0,
        no_numeric_validation_requested = TRUE
      )
    }

  }, error = function(e) {
    log_message("WARN", paste("Numeric data validation error:", e$message), category = "Validation")
    validation_result$numeric_data_summary <- list(
      validation_failed = TRUE,
      error_message = e$message
    )
  })

  # Step 3: Outlier assessment with safe operations
  log_message("DEBUG", "Outlier assessment", category = "Validation")

  tryCatch({
    outlier_summary <- list(
      total_outliers = 0,
      outlier_columns = 0
    )

    if (length(numeric_columns) > 0) {
      existing_numeric <- intersect(numeric_columns, names(data))

      for (col in existing_numeric) {
        col_data <- data[[col]]
        finite_data <- col_data[is.finite(col_data)]

        if (length(finite_data) > 4) {  # Need at least 5 points for IQR
          # Use simple IQR method for outlier detection
          Q1 <- quantile(finite_data, 0.25, na.rm = TRUE)
          Q3 <- quantile(finite_data, 0.75, na.rm = TRUE)
          IQR_val <- Q3 - Q1

          outliers <- finite_data < (Q1 - 1.5 * IQR_val) | finite_data > (Q3 + 1.5 * IQR_val)
          n_outliers <- sum(outliers, na.rm = TRUE)

          outlier_summary$total_outliers <- outlier_summary$total_outliers + n_outliers
          if (n_outliers > 0) {
            outlier_summary$outlier_columns <- outlier_summary$outlier_columns + 1
          }
        }
      }
    }

    validation_result$outlier_summary <- outlier_summary

  }, error = function(e) {
    log_message("WARN", paste("Outlier assessment error:", e$message), category = "Validation")
    validation_result$outlier_summary <- list(
      assessment_failed = TRUE,
      error_message = e$message
    )
  })

  # Step 4: Consistency checks
  log_message("DEBUG", "Consistency checks", category = "Validation")

  consistency_summary <- list(
    required_columns_present = TRUE,
    duplicate_rows = 0,
    consistency_score = 1.0
  )

  # Check required columns
  if (length(required_columns) > 0) {
    missing_required <- setdiff(required_columns, names(data))
    if (length(missing_required) > 0) {
      consistency_summary$required_columns_present <- FALSE
      consistency_summary$missing_required_columns <- missing_required
      consistency_summary$consistency_score <- 0.5
    }
  }

  # Check for duplicate rows (with error handling)
  tryCatch({
    consistency_summary$duplicate_rows <- sum(duplicated(data))
  }, error = function(e) {
    consistency_summary$duplicate_rows <- 0
  })

  validation_result$consistency_summary <- consistency_summary

  # Step 5: Calculate overall quality score with FIXED safe calculations
  tryCatch({
    # Component scores (0-1 scale)
    missing_score <- 1 - (validation_result$missing_data_summary$avg_missing_percentage %||% 0)

    numeric_score <- if (!is.null(validation_result$numeric_data_summary$avg_completeness)) {
      validation_result$numeric_data_summary$avg_completeness
    } else {
      0.8  # Default if no numeric data
    }

    consistency_score <- consistency_summary$consistency_score

    # Weighted overall score
    overall_score <- (missing_score * 0.4) + (numeric_score * 0.4) + (consistency_score * 0.2)

    # Determine quality grade
    quality_grade <- if (overall_score >= 0.9) {
      "Excellent"
    } else if (overall_score >= 0.8) {
      "Good"
    } else if (overall_score >= 0.7) {
      "Acceptable"
    } else if (overall_score >= 0.6) {
      "Needs Improvement"
    } else {
      "Poor"
    }

    validation_result$overall_quality <- list(
      score = overall_score,
      grade = quality_grade,
      validation_passed = overall_score >= (quality_thresholds$min_completeness %||% 0.7),
      component_scores = list(
        missing_data = missing_score,
        numeric_data = numeric_score,
        consistency = consistency_score
      )
    )

    log_message("INFO", paste("Overall data quality score:", round(overall_score, 2)), category = "Validation")
    log_message("INFO", paste("Quality grade:", quality_grade), category = "Validation")

  }, error = function(e) {
    log_message("WARN", paste("Quality score calculation error:", e$message), category = "Validation")

    # Fallback overall quality
    validation_result$overall_quality <- list(
      score = 0.8,  # Default decent score
      grade = "Good",
      validation_passed = TRUE,
      calculation_failed = TRUE,
      error_message = e$message
    )
  })

  return(validation_result)
}

#' Check Required Columns
#'
#' Validates that required columns exist and have appropriate data types.
#'
#' @param data Input data
#' @param column_specifications List of column specifications
#' @param strict_mode Whether to enforce strict type checking
#' @return Column validation results
#' @export
check_required_columns <- function(data, column_specifications, strict_mode = TRUE) {

  validation_results <- list(
    missing_columns = character(0),
    type_mismatches = character(0),
    validation_passed = TRUE
  )

  for (col_name in names(column_specifications)) {
    col_spec <- column_specifications[[col_name]]

    # Check existence
    if (!col_name %in% names(data)) {
      validation_results$missing_columns <- c(validation_results$missing_columns, col_name)
      validation_results$validation_passed <- FALSE
      next
    }

    # Check type if specified
    if (!is.null(col_spec$type) && strict_mode) {
      actual_type <- class(data[[col_name]])[1]
      expected_type <- col_spec$type

      if (actual_type != expected_type) {
        validation_results$type_mismatches <- c(
          validation_results$type_mismatches,
          paste0(col_name, " (expected: ", expected_type, ", actual: ", actual_type, ")")
        )
        validation_results$validation_passed <- FALSE
      }
    }
  }

  return(validation_results)
}

#' Validate Numeric Ranges
#'
#' Validates that numeric properties fall within realistic ranges.
#'
#' @param data Input data
#' @param property_ranges Named list of property ranges
#' @param action Action for out-of-range values: "warn", "error", or "clip"
#' @return Range validation results
#' @export
validate_numeric_ranges <- function(data, property_ranges, action = "warn") {

  validation_results <- list(
    range_violations = list(),
    corrected_data = data
  )

  for (prop in names(property_ranges)) {
    if (!prop %in% names(data)) {
      next
    }

    prop_data <- data[[prop]]
    range_spec <- property_ranges[[prop]]

    if (is.numeric(prop_data)) {
      violations <- which(prop_data < range_spec$min | prop_data > range_spec$max)

      if (length(violations) > 0) {
        validation_results$range_violations[[prop]] <- list(
          n_violations = length(violations),
          violation_indices = violations,
          out_of_range_values = prop_data[violations]
        )

        # Handle violations based on action
        if (action == "error") {
          stop("Range violations found in ", prop, ": ", length(violations), " values out of range")
        } else if (action == "clip") {
          validation_results$corrected_data[[prop]][violations] <- pmax(
            pmin(prop_data[violations], range_spec$max),
            range_spec$min
          )
        } else if (action == "warn") {
          log_message("WARN", paste("Range violations found in", prop, ":", length(violations), "values out of range"), category = "Validation")
        }
      }
    }
  }

  return(validation_results)
}

# ============================================================================
# 2. DATA TRANSFORMATION & PROCESSING UTILITIES
# ============================================================================

#' Standardize Property Names
#'
#' Standardizes soil property names across different data sources using
#' a comprehensive mapping system.
#'
#' @param data Input data with potentially inconsistent property names
#' @param property_mapping Custom property mapping (optional)
#' @param target_standard Target naming standard ("ssurgo", "nrcs", "custom")
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's log level so \code{INFO}-level progress messages print for the duration of this call (default \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @return Data with standardized property names
#' @export
standardize_property_names <- function(data,
                                       property_mapping = NULL,
                                       target_standard = "ssurgo",
                                       verbose = getOption("ssurgo.verbose", FALSE)) {
  .old_log_cfg <- set_verbose_logging(verbose)
  on.exit(options(soil_workflow_log_config = .old_log_cfg), add = TRUE)

  if (is.null(property_mapping)) {
    property_mapping <- get_default_property_mapping(target_standard)
  }

  log_message("INFO", paste("Standardizing property names to", target_standard, "standard"), category = "DataTransform")

  standardized_data <- data
  name_changes <- character(0)

  for (standard_name in names(property_mapping)) {
    source_names <- property_mapping[[standard_name]]

    # Find which source names exist in the data
    existing_sources <- intersect(source_names, names(data))

    if (length(existing_sources) > 0) {
      # Use the first available source (highest priority)
      source_name <- existing_sources[1]

      if (source_name != standard_name) {
        # Rename the column
        names(standardized_data)[names(standardized_data) == source_name] <- standard_name
        name_changes <- c(name_changes, paste(source_name, "->", standard_name))

        # Remove any additional source columns to avoid duplicates
        if (length(existing_sources) > 1) {
          additional_sources <- existing_sources[-1]
          standardized_data <- standardized_data[, !names(standardized_data) %in% additional_sources]
          name_changes <- c(name_changes, paste("Removed duplicates:", paste(additional_sources, collapse = ", ")))
        }
      }
    }
  }

  if (length(name_changes) > 0) {
    log_message("DEBUG", "Property name changes applied", category = "DataTransform")
    for (change in name_changes) {
      log_message("DEBUG", paste("  ", change), category = "DataTransform")
    }
  } else {
    log_message("DEBUG", "No property name changes needed", category = "DataTransform")
  }

  return(standardized_data)
}

#' Convert Depth Units
#'
#' Converts depth measurements between different units (cm, m, in, ft).
#'
#' @param depths Numeric vector of depth values
#' @param from_unit Source unit ("cm", "m", "in", "ft")
#' @param to_unit Target unit ("cm", "m", "in", "ft")
#' @param round_digits Number of decimal places for rounding
#' @return Converted depth values
#' @export
convert_depth_units <- function(depths, from_unit, to_unit, round_digits = 2) {

  if (from_unit == to_unit) {
    return(depths)
  }

  # Conversion factors to centimeters
  to_cm_factors <- list(
    "cm" = 1,
    "m" = 100,
    "in" = 2.54,
    "ft" = 30.48
  )

  if (!from_unit %in% names(to_cm_factors) || !to_unit %in% names(to_cm_factors)) {
    stop("Unsupported unit. Supported units: ", paste(names(to_cm_factors), collapse = ", "))
  }

  # Convert to cm first, then to target unit
  depths_cm <- depths * to_cm_factors[[from_unit]]
  converted_depths <- depths_cm / to_cm_factors[[to_unit]]

  # Round if specified
  if (!is.null(round_digits)) {
    converted_depths <- round(converted_depths, round_digits)
  }

  log_message("DEBUG", paste("Converted", length(depths), "depth values from", from_unit, "to", to_unit), category = "DataTransform")

  return(converted_depths)
}

#' Safe Coalesce
#'
#' Safely coalesces values from multiple columns, handling different data types
#' and missing values.
#'
#' @param data Input data
#' @param column_names Vector of column names in priority order
#' @param target_type Target data type ("numeric", "character", "logical")
#' @param default_value Default value if all columns are NA
#' @return Coalesced values
#' @export
safe_coalesce <- function(data, column_names, target_type = "numeric", default_value = NA) {

  # Check which columns exist
  existing_columns <- intersect(column_names, names(data))

  if (length(existing_columns) == 0) {
    log_message("WARN", paste("None of the specified columns exist:", paste(column_names, collapse = ", ")), category = "DataTransform")
    return(rep(default_value, nrow(data)))
  }

  result <- rep(default_value, nrow(data))

  # Coalesce in priority order
  for (col in existing_columns) {
    col_values <- data[[col]]

    # Type conversion if needed
    if (target_type == "numeric" && !is.numeric(col_values)) {
      col_values <- suppressWarnings(as.numeric(col_values))
    } else if (target_type == "character" && !is.character(col_values)) {
      col_values <- as.character(col_values)
    } else if (target_type == "logical" && !is.logical(col_values)) {
      col_values <- as.logical(col_values)
    }

    # Fill in missing values in result
    if (target_type == "numeric") {
      missing_in_result <- is.na(result)
      available_in_col <- !is.na(col_values)
      result[missing_in_result & available_in_col] <- col_values[missing_in_result & available_in_col]
    } else {
      missing_in_result <- is.na(result) | result == "" | result == default_value
      available_in_col <- !is.na(col_values) & col_values != ""
      result[missing_in_result & available_in_col] <- col_values[missing_in_result & available_in_col]
    }
  }

  n_filled <- sum(!is.na(result) & result != default_value)
  log_message("DEBUG", paste("Coalesced", n_filled, "values from columns:", paste(existing_columns, collapse = ", ")), category = "DataTransform")

  return(result)
}

#' Handle Missing Values
#'
#' Comprehensive missing value handling with multiple strategies.
#'
#' @param data Input data
#' @param strategy Missing value strategy ("remove", "interpolate", "mean", "median", "forward_fill")
#' @param columns Columns to process (NULL = all numeric columns)
#' @param group_by Grouping columns for grouped imputation
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's log level so \code{INFO}-level progress messages print for the duration of this call (default \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @return Data with handled missing values
#' @export
handle_missing_values <- function(data,
                                  strategy = "interpolate",
                                  columns = NULL,
                                  group_by = NULL,
                                  verbose = getOption("ssurgo.verbose", FALSE)) {
  .old_log_cfg <- set_verbose_logging(verbose)
  on.exit(options(soil_workflow_log_config = .old_log_cfg), add = TRUE)

  if (is.null(columns)) {
    columns <- names(data)[sapply(data, is.numeric)]
  }

  log_message("INFO", paste("Handling missing values using strategy:", strategy), category = "DataTransform")
  log_message("DEBUG", paste("Processing columns:", paste(columns, collapse = ", ")), category = "DataTransform")

  processed_data <- data

  for (col in columns) {
    if (!col %in% names(data)) {
      next
    }

    original_missing <- sum(is.na(data[[col]]))

    if (original_missing == 0) {
      next
    }

    if (!is.null(group_by) && all(group_by %in% names(data))) {
      # Grouped missing value handling
      processed_data <- processed_data |>
        group_by(across(all_of(group_by))) |>
        mutate(!!col := handle_missing_values_vector(.data[[col]], strategy)) |>
        ungroup()
    } else {
      # Global missing value handling
      processed_data[[col]] <- handle_missing_values_vector(data[[col]], strategy)
    }

    final_missing <- sum(is.na(processed_data[[col]]))
    filled <- original_missing - final_missing

    log_message("DEBUG", paste("  ", col, ": filled", filled, "of", original_missing, "missing values"), category = "DataTransform")
  }

  return(processed_data)
}

# ============================================================================
# 3. FILE I/O & DATA MANAGEMENT FUNCTIONS
# ============================================================================

#' Read Soil Data
#'
#' Robust soil data reading with automatic format detection and validation.
#'
#' @param file_path Path to data file
#' @param file_type File type ("auto", "csv", "xlsx", "rds", "txt")
#' @param validate_on_load Whether to validate data after loading
#' @param sheet_name Sheet name for Excel files
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's log level so \code{INFO}-level progress messages print for the duration of this call (default \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @return Loaded soil data
#' @export
read_soil_data <- function(file_path,
                           file_type = "auto",
                           validate_on_load = TRUE,
                           sheet_name = NULL,
                           verbose = getOption("ssurgo.verbose", FALSE)) {
  .old_log_cfg <- set_verbose_logging(verbose)
  on.exit(options(soil_workflow_log_config = .old_log_cfg), add = TRUE)

  if (!file.exists(file_path)) {
    stop("File not found: ", file_path)
  }

  log_message("INFO", paste("Reading soil data from:", file_path), category = "IO")

  # Auto-detect file type
  if (file_type == "auto") {
    file_ext <- tools::file_ext(file_path)
    file_type <- switch(tolower(file_ext),
                        "csv" = "csv",
                        "xlsx" = "xlsx",
                        "xls" = "xlsx",
                        "rds" = "rds",
                        "txt" = "txt",
                        "csv"  # default
    )
  }

  log_message("DEBUG", paste("Detected file type:", file_type), category = "IO")

  # Read data based on type
  data <- switch(file_type,
                 "csv" = read_csv_robust(file_path),
                 "xlsx" = read_excel_robust(file_path, sheet_name),
                 "rds" = readRDS(file_path),
                 "txt" = read_tsv_robust(file_path),
                 stop("Unsupported file type: ", file_type)
  )

  log_message("INFO", paste("Loaded data:", nrow(data), "rows,", ncol(data), "columns"), category = "IO")

  # Validate if requested
  if (validate_on_load) {
    log_message("DEBUG", "Validating loaded data", category = "IO")
    validation_results <- validate_data_quality(data)

    if (validation_results$overall_quality$score < 0.7) {
      log_message("WARN", paste("Loaded data has quality issues. Score:",
                                round(validation_results$overall_quality$score, 3)), category = "IO")
    }
  }

  return(data)
}

#' Write Soil Data
#'
#' Standardized soil data writing with metadata and backup options.
#'
#' @param data Data to write
#' @param file_path Output file path
#' @param file_type Output file type ("csv", "xlsx", "rds")
#' @param include_metadata Whether to include metadata
#' @param create_backup Whether to create backup
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's log level so \code{INFO}-level progress messages print for the duration of this call (default \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @return Success status
#' @export
write_soil_data <- function(data,
                            file_path,
                            file_type = "csv",
                            include_metadata = TRUE,
                            create_backup = TRUE,
                            verbose = getOption("ssurgo.verbose", FALSE)) {
  .old_log_cfg <- set_verbose_logging(verbose)
  on.exit(options(soil_workflow_log_config = .old_log_cfg), add = TRUE)

  log_message("INFO", paste("Writing soil data to:", file_path), category = "IO")
  log_message("DEBUG", paste("File type:", file_type), category = "IO")

  # Create backup if requested
  if (create_backup && file.exists(file_path)) {
    backup_path <- create_backup_filename(file_path)
    file.copy(file_path, backup_path)
    log_message("DEBUG", paste("Created backup:", backup_path), category = "IO")
  }

  # Create directory if needed
  dir.create(dirname(file_path), recursive = TRUE, showWarnings = FALSE)

  # Write data
  tryCatch({
    switch(file_type,
           "csv" = readr::write_csv(data, file_path),
           "xlsx" = write_excel_with_metadata(data, file_path, include_metadata),
           "rds" = saveRDS(data, file_path),
           stop("Unsupported file type: ", file_type)
    )

    log_message("INFO", "Data written successfully", category = "IO")
    return(TRUE)

  }, error = function(e) {
    log_message("ERROR", paste("Error writing data:", e$message), category = "IO")
    return(FALSE)
  })
}

#' Backup Data
#'
#' Creates versioned backups of important data files.
#'
#' @param source_path Source file path
#' @param backup_dir Backup directory (NULL = same directory)
#' @param max_backups Maximum number of backups to keep
#' @return Backup file path
#' @export
backup_data <- function(source_path, backup_dir = NULL, max_backups = 5) {

  if (!file.exists(source_path)) {
    stop("Source file not found: ", source_path)
  }

  if (is.null(backup_dir)) {
    backup_dir <- dirname(source_path)
  }

  dir.create(backup_dir, recursive = TRUE, showWarnings = FALSE)

  # Create timestamped backup filename
  base_name <- tools::file_path_sans_ext(basename(source_path))
  file_ext <- tools::file_ext(source_path)
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

  backup_filename <- paste0(base_name, "_backup_", timestamp, ".", file_ext)
  backup_path <- file.path(backup_dir, backup_filename)

  # Copy file
  if (file.copy(source_path, backup_path)) {
    log_message("DEBUG", paste("Backup created:", backup_path), category = "IO")

    # Clean up old backups
    cleanup_old_backups(backup_dir, base_name, max_backups)

    return(backup_path)
  } else {
    stop("Failed to create backup")
  }
}

#' Load Configuration
#'
#' Loads configuration from JSON or R files with validation and defaults.
#'
#' @param config_path Path to configuration file
#' @param config_type Configuration type ("json", "r", "auto")
#' @param validate_config Whether to validate configuration
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's log level so \code{INFO}-level progress messages print for the duration of this call (default \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @return Loaded configuration
#' @export
load_configuration <- function(config_path, config_type = "auto", validate_config = TRUE,
                               verbose = getOption("ssurgo.verbose", FALSE)) {
  .old_log_cfg <- set_verbose_logging(verbose)
  on.exit(options(soil_workflow_log_config = .old_log_cfg), add = TRUE)

  if (!file.exists(config_path)) {
    log_message("WARN", paste("Configuration file not found:", config_path, ". Using defaults."), category = "Config")
    return(get_default_configuration())
  }

  log_message("INFO", paste("Loading configuration from:", config_path), category = "Config")

  # Auto-detect config type
  if (config_type == "auto") {
    file_ext <- tools::file_ext(config_path)
    config_type <- switch(tolower(file_ext),
                          "json" = "json",
                          "r" = "r",
                          "json"  # default
    )
  }

  # Load configuration
  config <- switch(config_type,
                   "json" = jsonlite::fromJSON(config_path, simplifyVector = FALSE),
                   "r" = source(config_path)$value,
                   stop("Unsupported config type: ", config_type)
  )

  # Merge with defaults
  default_config <- get_default_configuration()
  config <- merge_configurations(default_config, config)

  # Validate if requested
  if (validate_config) {
    validation_results <- validate_configuration(config)
    if (!validation_results$valid) {
      log_message("WARN", paste("Configuration validation failed:", paste(validation_results$errors, collapse = ", ")), category = "Config")
    }
  }

  log_message("INFO", "Configuration loaded successfully", category = "Config")
  return(config)
}

# ============================================================================
# 4. STATISTICAL & MATHEMATICAL UTILITIES
# ============================================================================

#' Safe Correlation
#'
#' Numerically stable correlation calculation with handling for edge cases.
#'
#' @param x First variable
#' @param y Second variable (NULL for correlation matrix)
#' @param method Correlation method ("pearson", "spearman", "kendall")
#' @param handle_constant How to handle constant variables ("remove", "zero", "warn")
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's log level so \code{INFO}-level progress messages print for the duration of this call (default \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @return Correlation coefficient or matrix
#' @export
safe_correlation <- function(x, y = NULL, method = "pearson", handle_constant = "warn",
                             verbose = getOption("ssurgo.verbose", FALSE)) {
  .old_log_cfg <- set_verbose_logging(verbose)
  on.exit(options(soil_workflow_log_config = .old_log_cfg), add = TRUE)

  if (is.null(y)) {
    # Correlation matrix
    if (!is.data.frame(x) && !is.matrix(x)) {
      stop("For correlation matrix, x must be a data frame or matrix")
    }

    # Remove non-numeric columns
    if (is.data.frame(x)) {
      numeric_cols <- sapply(x, is.numeric)
      x <- x[, numeric_cols, drop = FALSE]
    }

    # Check for constant columns
    constant_cols <- apply(x, 2, function(col) var(col, na.rm = TRUE) == 0 | all(is.na(col)))

    if (any(constant_cols)) {
      constant_col_names <- names(x)[constant_cols]

      if (handle_constant == "remove") {
        x <- x[, !constant_cols, drop = FALSE]
        log_message("DEBUG", paste("Removed constant columns:", paste(constant_col_names, collapse = ", ")), category = "Statistics")
      } else if (handle_constant == "warn") {
        log_message("WARN", paste("Constant columns detected:", paste(constant_col_names, collapse = ", ")), category = "Statistics")
      }
    }

    if (ncol(x) < 2) {
      log_message("WARN", "Less than 2 columns remaining for correlation matrix", category = "Statistics")
      return(matrix(1, nrow = 1, ncol = 1, dimnames = list(names(x), names(x))))
    }

    # Calculate correlation matrix
    tryCatch({
      cor_matrix <- cor(x, method = method, use = "pairwise.complete.obs")

      # Check for numerical issues
      if (any(is.na(cor_matrix)) && method == "pearson") {
        log_message("WARN", "NAs in correlation matrix, trying spearman method", category = "Statistics")
        cor_matrix <- cor(x, method = "spearman", use = "pairwise.complete.obs")
      }

      return(cor_matrix)

    }, error = function(e) {
      log_message("WARN", paste("Correlation calculation failed:", e$message), category = "Statistics")
      n_vars <- ncol(x)
      return(diag(n_vars))
    })

  } else {
    # Pairwise correlation
    if (length(x) != length(y)) {
      stop("x and y must have the same length")
    }

    # Remove missing values
    complete_cases <- complete.cases(x, y)
    x_clean <- x[complete_cases]
    y_clean <- y[complete_cases]

    if (length(x_clean) < 3) {
      log_message("WARN", "Insufficient complete cases for correlation", category = "Statistics")
      return(NA)
    }

    # Check for constant variables
    if (var(x_clean) == 0 || var(y_clean) == 0) {
      if (handle_constant == "zero") {
        return(0)
      } else if (handle_constant == "warn") {
        log_message("WARN", "One or both variables are constant", category = "Statistics")
        return(0)
      } else {
        return(NA)
      }
    }

    # Calculate correlation
    tryCatch({
      cor_value <- cor(x_clean, y_clean, method = method)
      return(cor_value)

    }, error = function(e) {
      log_message("WARN", paste("Correlation calculation failed:", e$message), category = "Statistics")
      return(NA)
    })
  }
}

#' Calculate Confidence Intervals
#'
#' Calculates confidence intervals for various statistics.
#'
#' @param data Input data
#' @param statistic Statistic to calculate ("mean", "median", "proportion")
#' @param confidence_level Confidence level (default = 0.95)
#' @param method Method for calculation ("normal", "bootstrap", "t")
#' @param n_bootstrap Number of bootstrap samples
#' @return Confidence interval
#' @export
calculate_confidence_intervals <- function(data,
                                           statistic = "mean",
                                           confidence_level = 0.95,
                                           method = "normal",
                                           n_bootstrap = 1000) {

  if (any(is.na(data))) {
    data <- data[!is.na(data)]
    log_message("WARN", paste("Removed", sum(is.na(data)), "missing values"), category = "Statistics")
  }

  if (length(data) == 0) {
    return(list(lower = NA, upper = NA, estimate = NA))
  }

  alpha <- 1 - confidence_level

  if (method == "bootstrap") {
    # Bootstrap confidence intervals
    boot_stats <- replicate(n_bootstrap, {
      boot_sample <- sample(data, length(data), replace = TRUE)
      switch(statistic,
             "mean" = mean(boot_sample),
             "median" = median(boot_sample),
             "proportion" = mean(boot_sample),  # for binary data
             stop("Unsupported statistic: ", statistic)
      )
    })

    ci_lower <- quantile(boot_stats, alpha/2)
    ci_upper <- quantile(boot_stats, 1 - alpha/2)
    estimate <- switch(statistic,
                       "mean" = mean(data),
                       "median" = median(data),
                       "proportion" = mean(data)
    )

  } else if (method == "normal" || method == "t") {
    # Parametric confidence intervals
    n <- length(data)
    estimate <- mean(data)
    se <- sd(data) / sqrt(n)

    if (method == "t") {
      t_value <- qt(1 - alpha/2, df = n - 1)
      margin_error <- t_value * se
    } else {
      z_value <- qnorm(1 - alpha/2)
      margin_error <- z_value * se
    }

    ci_lower <- estimate - margin_error
    ci_upper <- estimate + margin_error

  } else {
    stop("Unsupported method: ", method)
  }

  return(list(
    lower = ci_lower,
    upper = ci_upper,
    estimate = estimate,
    confidence_level = confidence_level,
    method = method
  ))
}

#' Normalize Values
#'
#' Data normalization with multiple methods.
#'
#' @param x Input values
#' @param method Normalization method ("minmax", "zscore", "robust", "quantile")
#' @param center Center values (for zscore and robust)
#' @param scale Scale values (for zscore and robust)
#' @return Normalized values
#' @export
normalize_values <- function(x, method = "minmax", center = TRUE, scale = TRUE) {

  if (all(is.na(x))) {
    return(x)
  }

  switch(method,
         "minmax" = {
           min_x <- min(x, na.rm = TRUE)
           max_x <- max(x, na.rm = TRUE)
           if (min_x == max_x) {
             return(rep(0.5, length(x)))
           }
           (x - min_x) / (max_x - min_x)
         },

         "zscore" = {
           if (center && scale) {
             scale(x, center = TRUE, scale = TRUE)[, 1]
           } else if (center) {
             x - mean(x, na.rm = TRUE)
           } else if (scale) {
             x / sd(x, na.rm = TRUE)
           } else {
             x
           }
         },

         "robust" = {
           median_x <- median(x, na.rm = TRUE)
           mad_x <- mad(x, na.rm = TRUE)
           if (mad_x == 0) {
             return(x - median_x)
           }
           (x - median_x) / mad_x
         },

         "quantile" = {
           # Quantile normalization (rank-based)
           ranks <- rank(x, na.last = "keep")
           n_valid <- sum(!is.na(x))
           (ranks - 1) / (n_valid - 1)
         },

         stop("Unsupported normalization method: ", method)
  )
}

#' Detect Outliers
#'
#' Outlier detection using multiple methods.
#'
#' @param data Input data
#' @param method Detection method ("iqr", "zscore", "modified_zscore")
#' @param threshold Threshold parameter
#' @param return_indices Whether to return indices instead of logical vector
#' @return Outlier indicators or indices
#' @export
detect_outliers <- function(data, method = "iqr", threshold = NULL, return_indices = FALSE) {

  if (is.null(threshold)) {
    threshold <- switch(method,
                        "iqr" = 1.5,
                        "zscore" = 3,
                        "modified_zscore" = 3.5,
                        1.5
    )
  }

  outliers <- switch(method,
                     "iqr" = detect_outliers_iqr(data, threshold),
                     "zscore" = detect_outliers_zscore(data, threshold),
                     "modified_zscore" = detect_outliers_modified_zscore(data, threshold),
                     stop("Unsupported outlier detection method: ", method)
  )

  if (return_indices) {
    return(which(outliers))
  } else {
    return(outliers)
  }
}

# ============================================================================
# 5. LOGGING & ERROR HANDLING FUNCTIONS
# ============================================================================

#' Setup Logging
#'
#' Initializes comprehensive logging system for the workflow.
#'
#' @param log_file Log file path (NULL for console only)
#' @param log_level Log level ("DEBUG", "INFO", "WARN", "ERROR")
#' @param include_timestamp Whether to include timestamps
#' @param max_log_size Maximum log file size in MB
#' @return Logging configuration
#' @export
setup_logging <- function(log_file = NULL,
                          log_level = "INFO",
                          include_timestamp = TRUE,
                          max_log_size = 10) {

  # Create log directory if needed
  if (!is.null(log_file)) {
    dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)
  }

  # Initialize logging configuration
  log_config <- list(
    log_file = log_file,
    log_level = log_level,
    include_timestamp = include_timestamp,
    max_log_size = max_log_size,
    start_time = Sys.time()
  )

  # Set global logging configuration
  options(soil_workflow_log_config = log_config)

  # Log initialization
  log_message("INFO", "Logging system initialized", category = "System")
  log_message("INFO", paste("Log level:", log_level), category = "System")
  if (!is.null(log_file)) {
    log_message("INFO", paste("Log file:", log_file), category = "System")
  }

  return(log_config)
}

#' Log Message
#'
#' Comprehensive logging function with multiple levels and outputs.
#'
#' @param level Log level ("DEBUG", "INFO", "WARN", "ERROR")
#' @param message Log message
#' @param category Message category (optional)
#' @param details Additional details (optional)
#' @export
log_message <- function(level, message, category = NULL, details = NULL) {

  # Get logging configuration. Ambient fallback (used whenever setup_logging()
  # was never called) is deliberately "WARN", not "INFO" - the package is quiet
  # by default; INFO/DEBUG narration is opt-in via a function's own `verbose`
  # argument (see set_verbose_logging()) or an explicit setup_logging() call.
  log_config <- getOption("soil_workflow_log_config", list(
    log_file = NULL,
    log_level = "WARN",
    include_timestamp = TRUE
  ))

  # Check if message should be logged based on level
  level_hierarchy <- c("DEBUG" = 1, "INFO" = 2, "WARN" = 3, "ERROR" = 4)
  current_level <- level_hierarchy[log_config$log_level]
  message_level <- level_hierarchy[level]

  if (is.na(message_level) || message_level < current_level) {
    return(invisible())
  }

  # Format message
  if (log_config$include_timestamp) {
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    formatted_message <- paste0("[", timestamp, "] [", level, "] ", message)
  } else {
    formatted_message <- paste0("[", level, "] ", message)
  }

  # Add category if provided
  if (!is.null(category)) {
    formatted_message <- paste0(formatted_message, " [", category, "]")
  }

  # Add details if provided
  if (!is.null(details)) {
    formatted_message <- paste0(formatted_message, " - ", details)
  }

  # Output to console
  if (level %in% c("WARN", "ERROR")) {
    cat(formatted_message, "\n", file = stderr())
  } else {
    cat(formatted_message, "\n")
  }

  # Output to file if configured
  if (!is.null(log_config$log_file)) {
    # Check file size and rotate if needed
    if (file.exists(log_config$log_file)) {
      file_size_mb <- file.size(log_config$log_file) / (1024^2)
      if (file_size_mb > log_config$max_log_size) {
        rotate_log_file(log_config$log_file)
      }
    }

    # Write to log file
    cat(formatted_message, "\n", file = log_config$log_file, append = TRUE)
  }

  invisible()
}

#' Temporarily raise the package's log verbosity for the duration of a call
#'
#' Call at the top of a function with its own `verbose` argument, immediately followed by
#' `on.exit(options(soil_workflow_log_config = <the returned value>), add = TRUE)`. When
#' `verbose = TRUE`, raises the shared `soil_workflow_log_config` option's `log_level` to
#' `raised_level` for the remainder of the calling function's execution (restored via the
#' caller's own `on.exit()`), so every [log_message()] call made anywhere during that call -
#' including deep inside shared helpers that have no `verbose` parameter of their own, and
#' across package/file boundaries - becomes visible. Only ever *raises* the level, never
#' lowers it, so nested calls compose correctly regardless of order: an outer `verbose = TRUE`
#' makes everything below it visible even if an inner helper's own `verbose` is `FALSE`, and an
#' inner call's own restore is a harmless no-op when the level was already raised by an outer
#' caller. When `verbose = FALSE` (the default everywhere), does nothing - the ambient default
#' (see [log_message()]) is already quiet.
#'
#' @param verbose Logical.
#' @param raised_level Log level to switch to when `verbose = TRUE` (default `"INFO"` - matches
#'   the package's historical visible-by-default narration; `"DEBUG"`-level detail remains
#'   available only via an explicit `setup_logging(log_level = "DEBUG")` call, not through this
#'   per-call `verbose` mechanism).
#' @return Invisibly, the log config in effect *before* this call (pass to `on.exit()` to
#'   restore it).
#' @keywords internal
set_verbose_logging <- function(verbose, raised_level = "INFO") {
  # Fallback must match log_message()'s own fallback shape exactly (all 3 fields) -
  # on.exit() unconditionally writes this list back into the option even when
  # verbose = FALSE (restoring whatever was ambient before), so an incomplete list
  # here would silently replace a valid/unset option with a broken one, crashing
  # the next log_message() call on a missing field.
  old_config <- getOption("soil_workflow_log_config", list(
    log_file = NULL,
    log_level = "WARN",
    include_timestamp = TRUE
  ))
  if (isTRUE(verbose)) {
    options(soil_workflow_log_config = utils::modifyList(old_config, list(log_level = raised_level)))
  }
  invisible(old_config)
}

#' Handle Workflow Error
#'
#' Comprehensive error handling and recovery for workflow components.
#'
#' @param error Error object
#' @param context Context information
#' @param recovery_action Recovery action ("stop", "warn", "continue", "retry")
#' @param max_retries Maximum number of retries
#' @return Recovery action result
#' @export
handle_workflow_error <- function(error,
                                  context = "Unknown",
                                  recovery_action = "stop",
                                  max_retries = 3) {

  error_message <- paste("Error in", context, ":", error$message)
  log_message("ERROR", error_message, category = "ErrorHandler")

  # Log stack trace if available
  if (!is.null(error$call)) {
    log_message("DEBUG", paste("Call:", deparse(error$call)), category = "ErrorHandler")
  }

  # Handle based on recovery action
  switch(recovery_action,
         "stop" = {
           log_message("ERROR", "Stopping workflow due to error", category = "ErrorHandler")
           stop(error_message, call. = FALSE)
         },

         "warn" = {
           log_message("WARN", "Continuing with warning", category = "ErrorHandler")
           warning(error_message, call. = FALSE)
           return(list(status = "warning", message = error_message))
         },

         "continue" = {
           log_message("INFO", "Continuing despite error", category = "ErrorHandler")
           return(list(status = "error_ignored", message = error_message))
         },

         "retry" = {
           log_message("INFO", paste("Retrying operation, max retries:", max_retries), category = "ErrorHandler")
           return(list(status = "retry", message = error_message, max_retries = max_retries))
         },

         {
           log_message("ERROR", paste("Unknown recovery action:", recovery_action), category = "ErrorHandler")
           stop(error_message, call. = FALSE)
         }
  )
}

#' Track Progress
#'
#' Progress tracking utilities for long-running operations.
#'
#' @param current Current progress value
#' @param total Total expected value
#' @param message Progress message
#' @param update_frequency Update frequency (every nth call)
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's log level so \code{INFO}-level progress messages print for the duration of this call (default \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @export
track_progress <- function(current, total, message = "Processing", update_frequency = 10,
                           verbose = getOption("ssurgo.verbose", FALSE)) {
  .old_log_cfg <- set_verbose_logging(verbose)
  on.exit(options(soil_workflow_log_config = .old_log_cfg), add = TRUE)

  # Only update at specified frequency
  if (current %% update_frequency != 0 && current != total) {
    return(invisible())
  }

  percentage <- round((current / total) * 100, 1)
  progress_message <- paste0(message, ": ", current, "/", total, " (", percentage, "%)")

  # Estimate time remaining
  if (current > 0) {
    start_time <- getOption("progress_start_time", Sys.time())
    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    rate <- current / elapsed
    remaining <- (total - current) / rate

    if (remaining > 60) {
      time_str <- paste(round(remaining / 60, 1), "min remaining")
    } else {
      time_str <- paste(round(remaining, 1), "sec remaining")
    }

    progress_message <- paste(progress_message, "-", time_str)
  }

  log_message("INFO", progress_message, category = "Progress")

  # Set start time on first call
  if (current == 1) {
    options(progress_start_time = Sys.time())
  }

  invisible()
}

# ============================================================================
# 6. CONFIGURATION & PARAMETER MANAGEMENT
# ============================================================================

#' Get Default Configuration
#'
#' Returns comprehensive default configuration for the soil workflow.
#'
#' @param config_type Configuration type ("full", "minimal", "validation")
#' @return Default configuration
#' @export
get_default_configuration <- function(config_type = "full") {

  base_config <- list(
    workflow = list(
      version = "1.0.0",
      description = "Soil simulation workflow configuration"
    ),

    data_processing = list(
      max_depth = 250,
      exclude_unsuitable_horizons = TRUE,
      depth_units = "cm",
      missing_value_strategy = "interpolate"
    ),

    monte_carlo = list(
      n_simulations = 100,
      distribution_type = "triangular",
      confidence_level = 0.95,
      random_seed = 12345
    ),

    correlation = list(
      method = "pearson",
      minimum_samples = 30,
      handle_missing = "pairwise.complete.obs"
    ),

    gp_models = list(
      grouping_strategy = "auto",
      min_profiles_per_group = 3,
      min_observations_per_group = 15,
      optimize_hyperparameters = TRUE
    ),

    validation = list(
      run_validation = TRUE,
      quality_threshold = 0.7,
      generate_plots = TRUE,
      max_outlier_percentage = 5
    ),

    logging = list(
      log_level = "INFO",
      log_file = NULL,
      include_timestamp = TRUE
    )
  )

  if (config_type == "minimal") {
    # Return only essential configuration
    return(list(
      data_processing = base_config$data_processing,
      monte_carlo = list(n_simulations = base_config$monte_carlo$n_simulations),
      validation = list(run_validation = base_config$validation$run_validation)
    ))
  } else if (config_type == "validation") {
    # Return validation-focused configuration
    return(base_config$validation)
  } else {
    # Return full configuration
    return(base_config)
  }
}

#' Validate Parameters
#'
#' Validates input parameters for workflow functions.
#'
#' @param parameters Parameter list to validate
#' @param parameter_specs Parameter specifications
#' @param strict_mode Whether to use strict validation
#' @return Validation results
#' @export
validate_parameters <- function(parameters, parameter_specs, strict_mode = TRUE) {

  validation_result <- list(
    valid = TRUE,
    errors = character(0),
    warnings = character(0),
    corrected_parameters = parameters
  )

  # Validate each parameter
  for (param_name in names(parameter_specs)) {
    param_spec <- parameter_specs[[param_name]]
    param_value <- parameters[[param_name]]

    # Check if required parameter is missing
    if (is.null(param_value) || length(param_value) == 0) {
      if (param_spec$required %||% FALSE) {
        validation_result$errors <- c(validation_result$errors,
                                      paste("Required parameter missing:", param_name))
        validation_result$valid <- FALSE
      }
      next
    }

    # Type validation
    if (!is.null(param_spec$type)) {
      expected_type <- param_spec$type

      type_valid <- switch(expected_type,
                           "character" = is.character(param_value),
                           "numeric" = is.numeric(param_value),
                           "integer" = is.integer(param_value) || (is.numeric(param_value) && all(param_value == as.integer(param_value))),
                           "logical" = is.logical(param_value),
                           "list" = is.list(param_value),
                           TRUE  # Default to valid for unknown types
      )

      if (!type_valid) {
        validation_result$errors <- c(validation_result$errors,
                                      paste("Type mismatch for", param_name, "- expected", expected_type, "got", class(param_value)[1]))
        validation_result$valid <- FALSE
        next
      }
    }

    # Range validation (for numeric parameters)
    if (!is.null(param_spec$range) && is.numeric(param_value)) {
      range_limits <- param_spec$range

      # Check each value in the vector
      out_of_range <- param_value < range_limits[1] | param_value > range_limits[2]

      if (any(out_of_range, na.rm = TRUE)) {
        if (strict_mode) {
          validation_result$errors <- c(validation_result$errors,
                                        paste("Value out of range for", param_name, "- expected",
                                              paste(range_limits, collapse = " to ")))
          validation_result$valid <- FALSE
        } else {
          validation_result$warnings <- c(validation_result$warnings,
                                          paste("Value out of range for", param_name))
        }
      }
    }

    # Choices validation (FIXED to handle vectors properly)
    if (!is.null(param_spec$choices)) {
      valid_choices <- param_spec$choices

      # Handle vector parameters properly
      if (length(param_value) == 1) {
        # Single value - check directly
        if (!param_value %in% valid_choices) {
          validation_result$errors <- c(validation_result$errors,
                                        paste("Invalid choice for", param_name, "- got", param_value,
                                              "expected one of:", paste(valid_choices, collapse = ", ")))
          validation_result$valid <- FALSE
        }
      } else {
        # Vector values - check each element
        invalid_choices <- param_value[!param_value %in% valid_choices]

        if (length(invalid_choices) > 0) {
          validation_result$errors <- c(validation_result$errors,
                                        paste("Invalid choices for", param_name, "- got",
                                              paste(invalid_choices, collapse = ", "),
                                              "expected from:", paste(valid_choices, collapse = ", ")))
          validation_result$valid <- FALSE
        }
      }
    }

    # Length validation
    if (!is.null(param_spec$min_length)) {
      if (length(param_value) < param_spec$min_length) {
        validation_result$errors <- c(validation_result$errors,
                                      paste("Parameter", param_name, "too short - got", length(param_value),
                                            "expected at least", param_spec$min_length))
        validation_result$valid <- FALSE
      }
    }

    if (!is.null(param_spec$max_length)) {
      if (length(param_value) > param_spec$max_length) {
        if (strict_mode) {
          validation_result$errors <- c(validation_result$errors,
                                        paste("Parameter", param_name, "too long - got", length(param_value),
                                              "expected at most", param_spec$max_length))
          validation_result$valid <- FALSE
        } else {
          validation_result$warnings <- c(validation_result$warnings,
                                          paste("Parameter", param_name, "longer than recommended"))
        }
      }
    }
  }

  return(validation_result)
}

#' Merge Configurations
#'
#' Merges user configuration with default configuration.
#'
#' @param default_config Default configuration
#' @param user_config User-provided configuration
#' @param deep_merge Whether to perform deep merging
#' @return Merged configuration
#' @export
merge_configurations <- function(default_config, user_config, deep_merge = TRUE) {

  if (is.null(user_config)) {
    return(default_config)
  }

  if (!deep_merge) {
    # Simple merge - user config overrides defaults
    merged_config <- default_config
    for (name in names(user_config)) {
      merged_config[[name]] <- user_config[[name]]
    }
    return(merged_config)
  }

  # Deep merge - recursively merge nested lists
  merge_recursive <- function(default, user) {
    if (!is.list(default) || !is.list(user)) {
      return(user)  # User value overrides default
    }

    result <- default
    for (name in names(user)) {
      if (name %in% names(default) && is.list(default[[name]]) && is.list(user[[name]])) {
        result[[name]] <- merge_recursive(default[[name]], user[[name]])
      } else {
        result[[name]] <- user[[name]]
      }
    }
    return(result)
  }

  return(merge_recursive(default_config, user_config))
}

#' Export Workflow Metadata
#'
#' Exports comprehensive metadata about the workflow execution.
#'
#' @param workflow_results Workflow results
#' @param output_path Output file path
#' @param include_data_summary Whether to include data summaries
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's log level so \code{INFO}-level progress messages print for the duration of this call (default \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @return Success status
#' @export
export_workflow_metadata <- function(workflow_results,
                                     output_path,
                                     include_data_summary = TRUE,
                                     verbose = getOption("ssurgo.verbose", FALSE)) {
  .old_log_cfg <- set_verbose_logging(verbose)
  on.exit(options(soil_workflow_log_config = .old_log_cfg), add = TRUE)

  log_message("INFO", paste("Exporting workflow metadata to:", output_path), category = "Metadata")

  metadata <- list(
    workflow_info = list(
      timestamp = Sys.time(),
      r_version = R.version.string,
      platform = R.version$platform,
      packages = get_package_versions()
    ),

    execution_summary = extract_execution_summary(workflow_results),

    configuration = extract_configuration_info(workflow_results),

    quality_assessment = extract_quality_summary(workflow_results)
  )

  if (include_data_summary) {
    metadata$data_summary <- extract_data_summary(workflow_results)
  }

  # Export as JSON
  tryCatch({
    jsonlite::write_json(metadata, output_path, pretty = TRUE, auto_unbox = TRUE)
    log_message("INFO", "Metadata exported successfully", category = "Metadata")
    return(TRUE)

  }, error = function(e) {
    log_message("ERROR", paste("Error exporting metadata:", e$message), category = "Metadata")
    return(FALSE)
  })
}

# ============================================================================
# SUPPORTING HELPER FUNCTIONS
# ============================================================================

# Missing value handling for vectors
handle_missing_values_vector <- function(x, strategy) {
  switch(strategy,
         "remove" = x[!is.na(x)],
         "mean" = ifelse(is.na(x), mean(x, na.rm = TRUE), x),
         "median" = ifelse(is.na(x), median(x, na.rm = TRUE), x),
         "interpolate" = {
           if (all(is.na(x))) return(x)
           approx(seq_along(x), x, seq_along(x), method = "linear", rule = 2)$y
         },
         "forward_fill" = {
           for (i in 2:length(x)) {
             if (is.na(x[i]) && !is.na(x[i-1])) {
               x[i] <- x[i-1]
             }
           }
           x
         },
         x
  )
}

# File I/O helpers
read_csv_robust <- function(file_path) {
  tryCatch({
    readr::read_csv(file_path, show_col_types = FALSE)
  }, error = function(e) {
    # Fallback to base R
    read.csv(file_path, stringsAsFactors = FALSE)
  })
}

read_excel_robust <- function(file_path, sheet_name = NULL) {
  if (!requireNamespace("readxl", quietly = TRUE)) {
    stop("readxl package required for reading Excel files")
  }

  if (is.null(sheet_name)) {
    readxl::read_excel(file_path)
  } else {
    readxl::read_excel(file_path, sheet = sheet_name)
  }
}

read_tsv_robust <- function(file_path) {
  tryCatch({
    readr::read_tsv(file_path, show_col_types = FALSE)
  }, error = function(e) {
    read.delim(file_path, stringsAsFactors = FALSE)
  })
}

write_excel_with_metadata <- function(data, file_path, include_metadata) {
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    stop("openxlsx package required for writing Excel files")
  }

  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "Data")
  openxlsx::writeData(wb, "Data", data)

  if (include_metadata) {
    metadata <- list(
      "Export Date" = Sys.time(),
      "Rows" = nrow(data),
      "Columns" = ncol(data),
      "R Version" = R.version.string
    )

    openxlsx::addWorksheet(wb, "Metadata")
    openxlsx::writeData(wb, "Metadata",
                        data.frame(Field = names(metadata), Value = unlist(metadata)))
  }

  openxlsx::saveWorkbook(wb, file_path, overwrite = TRUE)
}

# Configuration helpers
get_default_property_mapping <- function(standard) {
  switch(standard,
         "ssurgo" = list(
           "claytotal_r" = c("clay_total", "clay_r", "clay", "claytotal"),
           "sandtotal_r" = c("sand_total", "sand_r", "sand", "sandtotal"),
           "silttotal_r" = c("silt_total", "silt_r", "silt", "silttotal"),
           "dbovendry_r" = c("bulk_density", "db_r", "db", "bulk_density_third_bar"),
           "ph1to1h2o_r" = c("ph", "pH", "ph_r", "pH_r"),
           "om_r" = c("om", "organic_matter", "organic_matter_r"),
           "cec7_r" = c("cec", "cec_r")
         ),
         "nrcs" = list(
           "clay_pct" = c("claytotal_r", "clay_total", "clay"),
           "sand_pct" = c("sandtotal_r", "sand_total", "sand"),
           "silt_pct" = c("silttotal_r", "silt_total", "silt"),
           "bulk_density" = c("dbovendry_r", "bulk_density_third_bar", "db"),
           "pH" = c("ph1to1h2o_r", "ph", "pH"),
           "organic_matter" = c("om_r", "om", "organic_matter"),
           "cec" = c("cec7_r", "cec")
         ),
         stop("Unsupported standard: ", standard)
  )
}

get_default_quality_thresholds <- function() {
  list(
    max_missing_percentage = 20,
    outlier_threshold = 3,
    min_sample_size = 10,
    correlation_threshold = 0.95
  )
}

# Data quality assessment helpers
assess_missing_data <- function(data, thresholds) {
  missing_summary <- data |>
    summarise_all(~ sum(is.na(.))) |>
    pivot_longer(everything(), names_to = "column", values_to = "missing_count") |>
    mutate(
      missing_percentage = .data$missing_count / nrow(data) * 100,
      quality_flag = .data$missing_percentage > thresholds$max_missing_percentage
    )

  high_missing_count <- sum(missing_summary$quality_flag)
  log_message("DEBUG", paste("Columns with >", thresholds$max_missing_percentage, "% missing:", high_missing_count), category = "DataAssessment")

  return(missing_summary)
}

validate_numeric_data <- function(data, numeric_columns, thresholds) {
  validation_results <- list()

  for (col in numeric_columns) {
    if (col %in% names(data)) {
      col_data <- data[[col]]

      validation_results[[col]] <- list(
        is_numeric = is.numeric(col_data),
        n_infinite = sum(is.infinite(col_data), na.rm = TRUE),
        n_negative = sum(col_data < 0, na.rm = TRUE),
        range = range(col_data, na.rm = TRUE)
      )
    }
  }

  return(validation_results)
}

assess_outliers <- function(data, thresholds) {
  numeric_cols <- names(data)[sapply(data, is.numeric)]
  outlier_summary <- list()

  for (col in numeric_cols) {
    outliers <- detect_outliers(data[[col]], method = "iqr")
    outlier_pct <- mean(outliers, na.rm = TRUE) * 100

    outlier_summary[[col]] <- list(
      n_outliers = sum(outliers, na.rm = TRUE),
      outlier_percentage = outlier_pct,
      quality_flag = outlier_pct > thresholds$max_missing_percentage
    )
  }

  return(outlier_summary)
}

perform_consistency_checks <- function(data, thresholds) {
  consistency_results <- list()

  # Check for duplicate rows
  n_duplicates <- sum(duplicated(data))
  consistency_results$duplicate_rows <- list(
    n_duplicates = n_duplicates,
    percentage = n_duplicates / nrow(data) * 100
  )

  # Check depth consistency
  if (all(c("hzdept_r", "hzdepb_r") %in% names(data))) {
    depth_issues <- sum(data$hzdept_r >= data$hzdepb_r, na.rm = TRUE)
    consistency_results$depth_consistency <- list(
      n_issues = depth_issues,
      percentage = depth_issues / nrow(data) * 100
    )
  }

  return(consistency_results)
}

calculate_overall_quality_score <- function(quality_results) {
  # Simple quality scoring based on missing data and outliers
  missing_score <- 1 - mean(quality_results$missing_data_assessment$missing_percentage, na.rm = TRUE) / 100

  if (!is.null(quality_results$outlier_assessment)) {
    outlier_scores <- sapply(quality_results$outlier_assessment, function(x) 1 - x$outlier_percentage / 100)
    outlier_score <- mean(outlier_scores, na.rm = TRUE)
  } else {
    outlier_score <- 1
  }

  overall_score <- mean(c(missing_score, outlier_score), na.rm = TRUE)

  grade <- if (overall_score >= 0.9) {
    "Excellent"
  } else if (overall_score >= 0.8) {
    "Good"
  } else if (overall_score >= 0.7) {
    "Acceptable"
  } else {
    "Poor"
  }

  return(list(score = overall_score, grade = grade))
}

# Outlier detection methods
detect_outliers_iqr <- function(x, threshold = 1.5) {
  q1 <- quantile(x, 0.25, na.rm = TRUE)
  q3 <- quantile(x, 0.75, na.rm = TRUE)
  iqr <- q3 - q1

  lower_bound <- q1 - threshold * iqr
  upper_bound <- q3 + threshold * iqr

  return(x < lower_bound | x > upper_bound)
}

detect_outliers_zscore <- function(x, threshold = 3) {
  z_scores <- abs((x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))
  return(z_scores > threshold)
}

detect_outliers_modified_zscore <- function(x, threshold = 3.5) {
  median_x <- median(x, na.rm = TRUE)
  mad_x <- mad(x, na.rm = TRUE)

  modified_z_scores <- 0.6745 * (x - median_x) / mad_x
  return(abs(modified_z_scores) > threshold)
}

# Backup and logging helpers
create_backup_filename <- function(file_path) {
  base_name <- tools::file_path_sans_ext(file_path)
  file_ext <- tools::file_ext(file_path)
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

  paste0(base_name, "_backup_", timestamp, ".", file_ext)
}

cleanup_old_backups <- function(backup_dir, base_name, max_backups) {
  # Find backup files for this base name
  backup_pattern <- paste0(base_name, "_backup_")
  backup_files <- list.files(backup_dir, pattern = backup_pattern, full.names = TRUE)

  if (length(backup_files) > max_backups) {
    # Sort by modification time and remove oldest
    backup_info <- file.info(backup_files)
    backup_info$filepath <- rownames(backup_info)
    backup_info <- backup_info[order(backup_info$mtime, decreasing = TRUE), ]

    files_to_remove <- backup_info$filepath[(max_backups + 1):nrow(backup_info)]

    for (file_path in files_to_remove) {
      if (file.remove(file_path)) {
        log_message("DEBUG", paste("Removed old backup:", basename(file_path)), category = "IO")
      }
    }
  }
}

rotate_log_file <- function(log_file) {
  backup_file <- paste0(tools::file_path_sans_ext(log_file), "_",
                        format(Sys.time(), "%Y%m%d_%H%M%S"), ".",
                        tools::file_ext(log_file))

  if (file.rename(log_file, backup_file)) {
    log_message("DEBUG", paste("Log file rotated to:", backup_file), category = "System")
  }
}

# Metadata extraction helpers
get_package_versions <- function() {
  installed_packages <- installed.packages()[, c("Package", "Version")]
  relevant_packages <- c("dplyr", "tidyr", "terra", "soilDB", "readr", "jsonlite")

  package_versions <- list()
  for (pkg in relevant_packages) {
    if (pkg %in% installed_packages[, "Package"]) {
      version <- installed_packages[installed_packages[, "Package"] == pkg, "Version"]
      package_versions[[pkg]] <- version
    }
  }

  return(package_versions)
}

extract_execution_summary <- function(workflow_results) {
  # Extract summary information from workflow results
  # Implementation depends on workflow results structure
  summary <- list(
    execution_time = "Not available",
    components_executed = "Not available",
    success_status = "Not available"
  )

  # Try to extract actual information if available
  if (!is.null(workflow_results) && is.list(workflow_results)) {
    if ("processing_metadata" %in% names(workflow_results)) {
      metadata <- workflow_results$processing_metadata
      if ("duration_seconds" %in% names(metadata)) {
        summary$execution_time <- paste(round(metadata$duration_seconds, 2), "seconds")
      }
      if ("input_rows" %in% names(metadata) && "output_rows" %in% names(metadata)) {
        summary$components_executed <- paste("Data processing:", metadata$input_rows, "->", metadata$output_rows, "rows")
      }
      summary$success_status <- "Completed"
    }
  }

  return(summary)
}

extract_configuration_info <- function(workflow_results) {
  # Extract configuration information
  config_info <- list(
    configuration_used = "Not available"
  )

  if (!is.null(workflow_results) && is.list(workflow_results)) {
    if ("processing_metadata" %in% names(workflow_results)) {
      metadata <- workflow_results$processing_metadata
      if ("options_used" %in% names(metadata)) {
        config_info$configuration_used <- metadata$options_used
      }
    }
  }

  return(config_info)
}

extract_quality_summary <- function(workflow_results) {
  # Extract quality assessment summary
  quality_summary <- list(
    overall_quality = "Not available"
  )

  if (!is.null(workflow_results) && is.list(workflow_results)) {
    if ("validation_results" %in% names(workflow_results)) {
      validation <- workflow_results$validation_results
      if (!is.null(validation) && "overall_quality" %in% names(validation)) {
        if ("score" %in% names(validation$overall_quality)) {
          quality_summary$overall_quality <- paste("Score:", round(validation$overall_quality$score, 3))
        }
      }
    }
    if ("quality_report" %in% names(workflow_results)) {
      quality <- workflow_results$quality_report
      if (!is.null(quality) && "overall_quality_score" %in% names(quality)) {
        quality_summary$overall_quality <- paste("Score:", round(quality$overall_quality_score, 3))
      }
    }
  }

  return(quality_summary)
}

extract_data_summary <- function(workflow_results) {
  # Extract data summary statistics
  data_summary <- list(
    data_dimensions = "Not available",
    property_summaries = "Not available"
  )

  if (!is.null(workflow_results) && is.list(workflow_results)) {
    if ("processed_data" %in% names(workflow_results)) {
      data <- workflow_results$processed_data
      if (is.data.frame(data)) {
        data_summary$data_dimensions <- paste(nrow(data), "rows x", ncol(data), "columns")

        # Summarize key soil properties if available
        soil_props <- names(data)[grepl("_r$", names(data))]
        if (length(soil_props) > 0) {
          prop_summary <- list()
          for (prop in soil_props[1:min(5, length(soil_props))]) {  # Limit to first 5
            prop_data <- data[[prop]]
            if (is.numeric(prop_data)) {
              prop_summary[[prop]] <- list(
                mean = round(mean(prop_data, na.rm = TRUE), 2),
                missing_pct = round(sum(is.na(prop_data)) / length(prop_data) * 100, 1)
              )
            }
          }
          data_summary$property_summaries <- prop_summary
        }
      }
    }
  }

  return(data_summary)
}

validate_configuration <- function(config) {
  # Basic configuration validation
  errors <- character(0)

  # Check required sections
  required_sections <- c("data_processing", "validation")
  missing_sections <- setdiff(required_sections, names(config))

  if (length(missing_sections) > 0) {
    errors <- c(errors, paste("Missing configuration sections:", paste(missing_sections, collapse = ", ")))
  }

  # Validate data processing section
  if ("data_processing" %in% names(config)) {
    dp_config <- config$data_processing

    if ("max_depth" %in% names(dp_config)) {
      if (!is.numeric(dp_config$max_depth) || dp_config$max_depth <= 0) {
        errors <- c(errors, "max_depth must be a positive number")
      }
    }

    if ("missing_value_strategy" %in% names(dp_config)) {
      valid_strategies <- c("remove", "interpolate", "mean", "median", "forward_fill")
      if (!dp_config$missing_value_strategy %in% valid_strategies) {
        errors <- c(errors, paste("Invalid missing_value_strategy. Valid options:", paste(valid_strategies, collapse = ", ")))
      }
    }
  }

  # Validate monte carlo section
  if ("monte_carlo" %in% names(config)) {
    mc_config <- config$monte_carlo

    if ("n_simulations" %in% names(mc_config)) {
      if (!is.numeric(mc_config$n_simulations) || mc_config$n_simulations < 1) {
        errors <- c(errors, "n_simulations must be a positive integer")
      }
    }

    if ("confidence_level" %in% names(mc_config)) {
      if (!is.numeric(mc_config$confidence_level) || mc_config$confidence_level <= 0 || mc_config$confidence_level >= 1) {
        errors <- c(errors, "confidence_level must be between 0 and 1")
      }
    }
  }

  # Validate logging section
  if ("logging" %in% names(config)) {
    log_config <- config$logging

    if ("log_level" %in% names(log_config)) {
      valid_levels <- c("DEBUG", "INFO", "WARN", "ERROR")
      if (!log_config$log_level %in% valid_levels) {
        errors <- c(errors, paste("Invalid log_level. Valid options:", paste(valid_levels, collapse = ", ")))
      }
    }
  }

  return(list(
    valid = length(errors) == 0,
    errors = errors
  ))
}

# ============================================================================
# SOIL PROPERTY VALIDATION
# ============================================================================

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
#' @export
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
  available_props <- get_available_properties(property_lookup)

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
get_available_properties <- function(property_lookup) {

  tryCatch({

    if (is.character(property_lookup) && length(property_lookup) == 1) {
      # Named source
      return(get_predefined_properties(property_lookup))

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
get_predefined_properties <- function(source_name) {

  switch(tolower(source_name),

         "ssurgo" = {
           # create_ssurgo_property_lookup_working() lives in ssurgo-acquisition.R
           # (migrated from mod01_ssurgo_data.R). The tryCatch() below is kept as
           # defensive degradation (WARN + empty vector) rather than removed,
           # since this switch statement has no other way to signal "unknown/
           # unavailable source" short of erroring.
           tryCatch({
             ssurgo_data <- create_ssurgo_property_lookup_working()
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
#' @export
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
#' Enhanced version that includes synonym matching
#'
#' @param properties Vector of property names to validate
#' @param property_lookup Property lookup source
#' @param strict_mode Enable strict validation
#' @param max_invalid_pct Maximum percentage of invalid properties allowed
#' @param performance_threshold Performance warning threshold
#'
#' @return Enhanced validation results with synonym matching
#'
#' @export
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

  # Define known SSURGO properties and their synonyms
  ssurgo_properties <- list(
    # Texture properties
    "sandtotal" = c("sandtotal", "sand_total", "sand", "sandtotal_r"),
    "claytotal" = c("claytotal", "clay_total", "clay", "claytotal_r"),
    "silttotal" = c("silttotal", "silt_total", "silt", "silttotal_r"),

    # Physical properties
    "dbovendry" = c("dbovendry", "bulk_density", "bd", "dbovendry_r"),
    "dbthirdbar" = c("dbthirdbar", "dbthirdbar_r"),

    # Chemical properties
    "ph1to1h2o" = c("ph1to1h2o", "ph", "ph_water", "ph1to1h2o_r"),
    "cec7" = c("cec7", "cec", "cation_exchange", "cec7_r"),
    "om" = c("om", "organic_matter", "som", "om_r"),

    # Water retention
    "wthirdbar" = c("wthirdbar", "water_33", "fc", "wthirdbar_r"),
    "wfifteenbar" = c("wfifteenbar", "water_1500", "pwp", "wfifteenbar_r"),
    "awc" = c("awc", "available_water", "awc_r"),

    # Other properties
    "rfv" = c("rfv", "rock_fragments", "fragvol", "rfv_r"),
    "ksat" = c("ksat", "sat_hydraulic_cond", "ksat_r"),
    "ec" = c("ec", "electrical_conductivity", "ec_r")
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
#' @export
get_default_synonyms <- function(property_lookup) {

  if (is.character(property_lookup) && length(property_lookup) == 1) {
    source_name <- tolower(property_lookup)
  } else {
    source_name <- "general"
  }

  switch(source_name,

         "ssurgo" = list(
           "claytotal_r" = c("clay", "clay_total", "clay_pct"),
           "sandtotal_r" = c("sand", "sand_total", "sand_pct"),
           "silttotal_r" = c("silt", "silt_total", "silt_pct"),
           "dbovendry_r" = c("bulk_density", "bd", "bulk_density_third_bar"),
           "ph1to1h2o_r" = c("ph", "pH", "ph_water"),
           "om_r" = c("om", "organic_matter"),
           "cec7_r" = c("cec", "cation_exchange_capacity")
         ),

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
#' @export
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
      # validate_properties_with_synonyms() has no `synonyms` parameter (it
      # uses its own internal SSURGO synonym table) - this call previously
      # passed lookup$synonyms positionally into strict_mode's slot and
      # strict_mode into max_invalid_pct's slot, silently mis-validating
      # whenever a property_lookup object had synonyms set. Pass strict_mode
      # by name to the actual signature instead.
      return(validate_properties_with_synonyms(props, lookup$properties, strict_mode = strict_mode))
    }
  }

  class(lookup) <- "property_lookup"
  return(lookup)
}

# ============================================================================
# GEOMETRY VALIDATION
# ============================================================================

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
#' @export
validate_wkt_geometry <- function(wkt_string,
                                  crs = "epsg:4326",
                                  validation_context = "geographic",
                                  area_limits = list(min = 0.0001, max = 100, max_action = "warn"),
                                  complexity_limits = list(max_vertices = 1000, max_parts = 100),
                                  bounds_check = list(xmin = -180, xmax = 180, ymin = -90, ymax = 90),
                                  strict_mode = TRUE,
                                  return_geometry = FALSE) {

  log_message("DEBUG", "Starting WKT validation in geographic context", category = "GeometryValidation")

  validation_result <- list(
    valid = TRUE,
    errors = character(0),
    warnings = character(0),
    geometry_stats = list()
  )

  tryCatch({
    log_message("DEBUG", paste("Parsing WKT geometry with CRS:", crs), category = "GeometryValidation")

    # Parse WKT with proper error handling
    if (!requireNamespace("terra", quietly = TRUE)) {
      stop("terra package required for geometry validation")
    }

    # Create geometry object
    geom <- terra::vect(wkt_string, crs = crs)

    # Check if geometry is valid
    if (!terra::is.valid(geom)) {
      validation_result$errors <- c(validation_result$errors, "Invalid WKT geometry")
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

  }, error = function(e) {
    log_message("ERROR", paste("WKT parsing error:", e$message), category = "GeometryValidation")
    validation_result$warnings <- c(validation_result$warnings, paste("Geometry validation error:", e$message))

    # Provide default geometry stats to prevent downstream errors
    validation_result$geometry_stats <- list(
      bbox = c(-180, -90, 180, 90),
      area_deg2 = 1.0,
      area_units = "degrees^2",
      geometry_type = "POLYGON",
      n_vertices = 4
    )
  })

  log_message("DEBUG", paste("WKT validation complete - Valid:", validation_result$valid,
                             "Errors:", length(validation_result$errors),
                             "Warnings:", length(validation_result$warnings)),
              category = "GeometryValidation")

  return(validation_result)
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
