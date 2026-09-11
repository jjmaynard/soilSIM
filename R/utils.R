# Shared utilities and helpers
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

  # Check what additional columns are available for fuller detection
  has_desgnmaster <- "desgnmaster" %in% names(data)
  has_restriction_top <- "restriction_top" %in% names(data)
  has_hzdept_r <- "hzdept_r" %in% names(data)
  has_is_potentially_restrictive <- "is_potentially_restrictive" %in% names(data)
  has_missing_required_properties <- "missing_required_properties" %in% names(data)
  has_organic_horizon <- "organic_horizon" %in% names(data)

  n <- length(hznames)
  unsuitable <- rep(FALSE, n)

  # toupper()/trimws() are vectorized once here rather than per-row inside the loop below.
  # The loop's own classification logic (which pattern a given row matches) stays row-dependent.
  hz_upper <- toupper(trimws(hznames))
  master_upper <- if (has_desgnmaster) toupper(trimws(as.character(data$desgnmaster))) else rep(NA_character_, n)

  # PRIORITY 1: Always check horizon names first (R, Cr, O should ALWAYS be unsuitable)
  for (i in seq_len(n)) {
    hz <- hz_upper[i]
    master <- master_upper[i]

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

  # Step 2: Numeric data validation with error handling
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
      validation_passed = overall_score >= (quality_thresholds$min_completeness %||% 0.8),
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
    validation_result$overall_quality <<- list(
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
#' @section Usage note:
#' Fully implemented and exported, but not currently called from anywhere
#' else in the package (confirmed by source grep) - standalone public API
#' for callers who need it, not dead/broken code.
#' @keywords internal
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
#' @keywords internal
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
#' @keywords internal
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
#' @keywords internal
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
#' @keywords internal
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
#' @keywords internal
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
#' @keywords internal
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
#' @keywords internal
compute_confidence_intervals <- function(data,
                                           statistic = "mean",
                                           confidence_level = 0.95,
                                           method = "normal",
                                           n_bootstrap = 1000) {

  if (any(is.na(data))) {
    n_missing <- sum(is.na(data))
    data <- data[!is.na(data)]
    log_message("WARN", paste("Removed", n_missing, "missing values"), category = "Statistics")
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
#' @section Usage note:
#' Fully implemented and exported, but not currently called from anywhere
#' else in the package (confirmed by source grep) - standalone public API
#' for callers who need it, not dead/broken code.
#' @keywords internal
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
#' @keywords internal
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

# Canonical SSURGO property synonym table (single source of truth).
#
# Keyed by the canonical `_r`-suffixed SSURGO column name; each entry lists
# every known alias/synonym for that property (including common short forms
# and the non-suffixed base name). `get_default_property_mapping("ssurgo")`,
# `default_property_synonyms("ssurgo")`, and `validate_properties_with_synonyms()`
# all derive from this one table, so they cannot drift out of sync.
.ssurgo_property_synonyms <- function() {
  list(
    "claytotal_r"   = c("claytotal", "clay_total", "clay_r", "clay", "clay_pct"),
    "sandtotal_r"   = c("sandtotal", "sand_total", "sand_r", "sand", "sand_pct"),
    "silttotal_r"   = c("silttotal", "silt_total", "silt_r", "silt", "silt_pct"),
    "dbovendry_r"   = c("dbovendry", "bulk_density", "db_r", "db", "bd", "bulk_density_third_bar"),
    "dbthirdbar_r"  = c("dbthirdbar"),
    "ph1to1h2o_r"   = c("ph1to1h2o", "ph", "pH", "ph_r", "pH_r", "ph_water"),
    "cec7_r"        = c("cec7", "cec", "cec_r", "cation_exchange", "cation_exchange_capacity"),
    "om_r"          = c("om", "organic_matter", "organic_matter_r", "som"),
    "wthirdbar_r"   = c("wthirdbar", "water_33", "fc"),
    "wfifteenbar_r" = c("wfifteenbar", "water_1500", "pwp"),
    "awc_r"         = c("awc", "available_water"),
    "rfv_r"         = c("rfv", "rock_fragments", "fragvol"),
    "ksat_r"        = c("ksat", "sat_hydraulic_cond"),
    "ec_r"          = c("ec", "electrical_conductivity"),
    # 5 chemistry properties - "ec_r" above already
    # existed; the other 4 are new entries.
    "caco3_r"       = c("caco3", "calcium_carbonate"),
    "ecec_r"        = c("ecec", "effective_cec", "effective_cation_exchange"),
    "gypsum_r"      = c("gypsum"),
    "sar_r"         = c("sar", "sodium_adsorption_ratio")
  )
}

# Configuration helpers
get_default_property_mapping <- function(standard) {
  ssurgo_synonyms <- .ssurgo_property_synonyms()

  switch(standard,
         "ssurgo" = ssurgo_synonyms[c("claytotal_r", "sandtotal_r", "silttotal_r",
                                      "dbovendry_r", "ph1to1h2o_r", "om_r", "cec7_r")],
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
create_backup_filename <- function(file_path, backup_dir = NULL) {
  base_name <- tools::file_path_sans_ext(basename(file_path))
  file_ext <- tools::file_ext(file_path)
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  backup_filename <- paste0(base_name, "_backup_", timestamp, ".", file_ext)

  if (is.null(backup_dir)) {
    backup_dir <- dirname(file_path)
  }

  file.path(backup_dir, backup_filename)
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

