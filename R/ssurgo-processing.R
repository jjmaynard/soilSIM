#' @title SSURGO Data Processing & Standardization
#' @description Functions for processing raw SSURGO data into standardized formats
#' @name ssurgo_processing
NULL

# ==============================================================================
# MAIN PROCESSING FUNCTIONS
# ==============================================================================

#' Process SSURGO Data (Main Entry Point)
#'
#' Enhanced version that maintains compatibility with working infill_soil_property()
#' while adding comprehensive processing, validation, and reporting capabilities.
#' Uses Module 8 utilities for general processing tasks.
#'
#' @param raw_data Raw SSURGO data from download functions
#' @param processing_options List of processing options and parameters
#' @param validate_results Logical; validate processing results (default: TRUE)
#' @param max_depth Maximum depth (cm) for processing (default: 250)
#' @param verbose Logical; provide detailed progress messages
#'
#' @return List containing:
#'   \itemize{
#'     \item processed_data: Combined processed data (compatible with infill functions)
#'     \item horizon_data: Processed horizon data
#'     \item component_data: Processed component data
#'     \item processing_metadata: Processing information and statistics
#'     \item validation_results: Data validation results
#'     \item quality_report: Comprehensive quality assessment
#'   }
#'
#' @export
process_ssurgo_data <- function(raw_data,
                                processing_options = list(),
                                validate_results = TRUE,
                                max_depth = 250,
                                verbose = getOption("ssurgo.verbose", FALSE)) {
  .old_log_cfg <- set_verbose_logging(verbose)
  on.exit(options(soil_workflow_log_config = .old_log_cfg), add = TRUE)

  start_time <- Sys.time()

  if (verbose) log_message("INFO", "SSURGO Data Processing Started", category = "Processing")

  # Set default processing options
  default_options <- list(
    detect_unsuitable = TRUE,         # Use Module 8 is_unsuitable function
    advanced_cleaning = TRUE,         # Use enhanced cleaning functions
    standardize_names = TRUE,         # Column name standardization using Module 8
    remove_invalid = TRUE,           # Remove invalid records
    calculate_derived = TRUE,        # Calculate derived properties
    validate_logic = TRUE,           # Data validation using Module 8
    preserve_compatibility = TRUE    # Maintain compatibility with infill functions
  )

  options <- merge_configurations(default_options, processing_options)

  if (verbose) {
    log_message("INFO", "Processing options configured", category = "Processing")
    log_message("DEBUG", paste("Input data:", nrow(raw_data), "rows"), category = "Processing")
  }

  # Step 1: Process horizon data using working-compatible logic
  if (verbose) log_message("INFO", "Processing horizon data", category = "Processing")

  horizon_result <- process_horizon_data_working_compatible(
    raw_data = raw_data,
    detect_unsuitable = options$detect_unsuitable,
    advanced_cleaning = options$advanced_cleaning,
    standardize_names = options$standardize_names,
    remove_invalid = options$remove_invalid,
    calculate_derived = options$calculate_derived,
    max_depth = max_depth,
    verbose = verbose
  )

  # Step 2: Process component data
  if (verbose) log_message("INFO", "Processing component data", category = "Processing")

  component_result <- process_component_data_working_compatible(
    raw_data = raw_data,
    standardize_names = options$standardize_names,
    remove_invalid = options$remove_invalid,
    verbose = verbose
  )

  # Step 3: Create main processed dataset (compatible with infill functions)
  if (verbose) log_message("INFO", "Creating compatible processed dataset", category = "Processing")

  processed_data <- create_infill_compatible_dataset(
    raw_data = raw_data,
    horizon_processing = horizon_result,
    component_processing = component_result,
    options = options,
    verbose = verbose
  )

  # Step 4: Validate processing results using Module 8
  validation_results <- NULL
  quality_report <- NULL

  if (validate_results || options$validate_logic) {
    if (verbose) log_message("INFO", "Validating processing results", category = "Processing")

    # Use Module 8 validation functions
    validation_results <- validate_data_quality(
      processed_data,
      required_columns = c("cokey", "hzname", "hzdept_r", "hzdepb_r"),
      numeric_columns = names(processed_data)[grepl("_r$", names(processed_data))]
    )

    # Generate comprehensive quality report
    quality_report <- generate_processing_quality_report(
      original_data = raw_data,
      processed_data = processed_data,
      horizon_stats = horizon_result$processing_stats,
      component_stats = component_result$processing_stats,
      validation_results = validation_results
    )
  }

  # Step 5: Generate processing metadata
  end_time <- Sys.time()
  processing_metadata <- list(
    processing_start = start_time,
    processing_end = end_time,
    duration_seconds = as.numeric(difftime(end_time, start_time, units = "secs")),
    options_used = options,
    max_depth = max_depth,
    input_rows = nrow(raw_data),
    output_rows = nrow(processed_data),
    unique_cokeys = length(unique(processed_data$cokey %||% character(0))),
    compatible_with_infill = options$preserve_compatibility,
    horizon_processing_stats = horizon_result$processing_stats,
    component_processing_stats = component_result$processing_stats
  )

  if (verbose) {
    log_message("INFO", "Processing Complete", category = "Processing")
    log_message("INFO", paste("Duration:", round(processing_metadata$duration_seconds, 2), "seconds"), category = "Processing")
    log_message("INFO", paste("Input rows:", processing_metadata$input_rows), category = "Processing")
    log_message("INFO", paste("Output rows:", processing_metadata$output_rows), category = "Processing")

    if ("unsuitable_horizon" %in% names(processed_data)) {
      unsuitable_count <- sum(processed_data$unsuitable_horizon, na.rm = TRUE)
      log_message("INFO", paste("Unsuitable horizons flagged:", unsuitable_count), category = "Processing")
    }
  }

  # Return comprehensive results with main processed_data first for compatibility
  return(list(
    processed_data = processed_data,           # Main output (compatible with infill functions)
    horizon_data = horizon_result$processed_data,
    component_data = component_result$processed_data,
    processing_metadata = processing_metadata,
    validation_results = validation_results,
    quality_report = quality_report
  ))
}

#' Process Horizon Data (Working Compatible)
#'
#' Enhanced horizon processing that uses Module 8 utilities and proven working functions
#'
#' @param raw_data Raw SSURGO data containing horizon information
#' @param detect_unsuitable Logical; detect unsuitable horizons using Module 8 function
#' @param advanced_cleaning Logical; apply advanced data cleaning functions
#' @param standardize_names Logical; standardize column names using Module 8
#' @param remove_invalid Logical; remove invalid horizon records
#' @param calculate_derived Logical; calculate derived properties
#' @param max_depth Maximum depth for processing
#' @param verbose Logical; provide progress messages
#'
#' @return List with processed horizon data and processing statistics
#'
#' @export
process_horizon_data_working_compatible <- function(raw_data,
                                                    detect_unsuitable = TRUE,
                                                    advanced_cleaning = TRUE,
                                                    standardize_names = TRUE,
                                                    remove_invalid = TRUE,
                                                    calculate_derived = TRUE,
                                                    max_depth = 250,
                                                    verbose = FALSE) {

  processing_stats <- list()
  cleaning_reports <- list()

  # Start with copy of raw data
  horizon_data <- raw_data
  processing_stats$initial_rows <- nrow(horizon_data)

  if (verbose) log_message("DEBUG", paste("Initial horizon records:", processing_stats$initial_rows), category = "HorizonProcessing")

  # Step 1: Standardize column names using Module 8
  if (standardize_names) {
    if (verbose) log_message("DEBUG", "Standardizing column names", category = "HorizonProcessing")
    horizon_data <- standardize_property_names(horizon_data, target_standard = "ssurgo")
  }

  # Step 2: Detect unsuitable horizons using Module 8 function
  if (detect_unsuitable) {
    if (verbose) log_message("DEBUG", "Detecting unsuitable horizons", category = "HorizonProcessing")
    horizon_data$unsuitable_horizon <- is_unsuitable(horizon_data, hzname_col = "hzname")

    unsuitable_count <- sum(horizon_data$unsuitable_horizon, na.rm = TRUE)
    processing_stats$unsuitable_horizons <- unsuitable_count

    if (verbose && unsuitable_count > 0) {
      log_message("DEBUG", paste("Found", unsuitable_count, "unsuitable horizons"), category = "HorizonProcessing")
    }
  }

  # Step 3: Advanced data cleaning using SSURGO-specific functions
  if (advanced_cleaning) {
    if (verbose) log_message("DEBUG", "Applying advanced data cleaning", category = "HorizonProcessing")

    # Identify soil property columns
    property_columns <- identify_soil_property_columns_working(horizon_data)

    for (prop_base in property_columns) {
      if (verbose) log_message("DEBUG", paste("Cleaning", prop_base), category = "HorizonProcessing")

      cleaning_result <- clean_property_data_ssurgo_compatible(
        horizon_data, prop_base,
        generate_report = TRUE,
        verbose = FALSE
      )

      horizon_data <- cleaning_result$data
      cleaning_reports[[prop_base]] <- cleaning_result$report
    }

    processing_stats$cleaning_reports <- cleaning_reports
  }

  # Step 4: Remove invalid horizon records
  if (remove_invalid) {
    if (verbose) log_message("DEBUG", "Removing invalid horizon records", category = "HorizonProcessing")

    initial_count <- nrow(horizon_data)
    horizon_data <- remove_invalid_horizons_working_compatible(horizon_data, max_depth = max_depth, verbose = verbose)
    processing_stats$removed_invalid <- initial_count - nrow(horizon_data)

    if (verbose && processing_stats$removed_invalid > 0) {
      log_message("DEBUG", paste("Removed", processing_stats$removed_invalid, "invalid records"), category = "HorizonProcessing")
    }
  }

  # Step 5: Calculate derived properties
  if (calculate_derived) {
    if (verbose) log_message("DEBUG", "Calculating derived horizon properties", category = "HorizonProcessing")
    horizon_data <- calculate_derived_horizon_properties_working(horizon_data, verbose = verbose)
  }

  # Step 6: Calculate processing statistics
  processing_stats$final_rows <- nrow(horizon_data)
  processing_stats$rows_retained <- if (processing_stats$initial_rows > 0) {
    processing_stats$final_rows / processing_stats$initial_rows
  } else {
    0
  }
  processing_stats$unique_cokeys <- length(unique(horizon_data$cokey %||% character(0)))

  if ("hzdept_r" %in% names(horizon_data) && "hzdepb_r" %in% names(horizon_data)) {
    processing_stats$depth_range <- range(c(horizon_data$hzdept_r, horizon_data$hzdepb_r), na.rm = TRUE)
    processing_stats$avg_horizon_thickness <- mean(horizon_data$hzdepb_r - horizon_data$hzdept_r, na.rm = TRUE)
  }

  # Calculate property completeness
  processing_stats$property_completeness <- calculate_property_completeness_working(horizon_data)

  if (verbose) {
    log_message("DEBUG", paste("Final horizon records:", processing_stats$final_rows), category = "HorizonProcessing")
    log_message("DEBUG", paste("Retention rate:", round(processing_stats$rows_retained * 100, 1), "%"), category = "HorizonProcessing")
  }

  return(list(
    processed_data = horizon_data,
    processing_stats = processing_stats
  ))
}

#' Process Component Data (Working Compatible)
#'
#' Enhanced component processing using Module 8 utilities
#'
#' @param raw_data Raw SSURGO data containing component information
#' @param standardize_names Logical; standardize column names using Module 8
#' @param remove_invalid Logical; remove invalid component records
#' @param verbose Logical; provide progress messages
#'
#' @return List with processed component data and processing statistics
#'
#' @export
process_component_data_working_compatible <- function(raw_data,
                                                      standardize_names = TRUE,
                                                      remove_invalid = TRUE,
                                                      verbose = FALSE) {

  processing_stats <- list()

  # Extract component-level information
  component_columns <- c("cokey", "compname", "comppct_r", "majcompflag",
                         "localphase", "slope_r", "taxclname", "taxorder",
                         "taxsuborder", "taxgrtgroup", "taxsubgrp", "taxpartsize",
                         "mukey", "compkind", "hydgrp", "hydric")

  available_columns <- intersect(component_columns, names(raw_data))

  if (length(available_columns) == 0) {
    log_message("WARN", "No component columns found in raw data", category = "ComponentProcessing")
    return(list(
      processed_data = data.frame(),
      processing_stats = list(initial_rows = 0, final_rows = 0)
    ))
  }

  # Create component data by selecting unique component records
  if ("cokey" %in% available_columns) {
    component_data <- raw_data |>
      dplyr::select(dplyr::all_of(available_columns)) |>
      dplyr::distinct()
  } else {
    component_data <- raw_data |>
      dplyr::select(dplyr::all_of(available_columns))
  }

  processing_stats$initial_rows <- nrow(component_data)

  if (verbose) log_message("DEBUG", paste("Initial component records:", processing_stats$initial_rows), category = "ComponentProcessing")

  # Step 1: Standardize column names using Module 8
  if (standardize_names) {
    if (verbose) log_message("DEBUG", "Standardizing component column names", category = "ComponentProcessing")
    component_data <- standardize_property_names(component_data, target_standard = "ssurgo")
  }

  # Step 2: Remove invalid component records
  if (remove_invalid) {
    if (verbose) log_message("DEBUG", "Removing invalid component records", category = "ComponentProcessing")

    initial_count <- nrow(component_data)
    component_data <- remove_invalid_components_working(component_data, verbose = verbose)
    processing_stats$removed_invalid <- initial_count - nrow(component_data)

    if (verbose && processing_stats$removed_invalid > 0) {
      log_message("DEBUG", paste("Removed", processing_stats$removed_invalid, "invalid records"), category = "ComponentProcessing")
    }
  }

  # Step 3: Calculate basic component statistics
  if (verbose) log_message("DEBUG", "Calculating component statistics", category = "ComponentProcessing")
  component_data <- calculate_component_stats_working(component_data, verbose = verbose)

  # Step 4: Calculate processing statistics
  processing_stats$final_rows <- nrow(component_data)
  processing_stats$rows_retained <- if (processing_stats$initial_rows > 0) {
    processing_stats$final_rows / processing_stats$initial_rows
  } else {
    0
  }
  processing_stats$unique_cokeys <- length(unique(component_data$cokey %||% character(0)))

  # Component-specific statistics
  if ("comppct_r" %in% names(component_data)) {
    processing_stats$component_percent_range <- range(component_data$comppct_r, na.rm = TRUE)
    processing_stats$avg_component_percent <- mean(component_data$comppct_r, na.rm = TRUE)
  }

  if (verbose) {
    log_message("DEBUG", paste("Final component records:", processing_stats$final_rows), category = "ComponentProcessing")
    log_message("DEBUG", paste("Retention rate:", round(processing_stats$rows_retained * 100, 1), "%"), category = "ComponentProcessing")
  }

  return(list(
    processed_data = component_data,
    processing_stats = processing_stats
  ))
}

#' Create Infill-Compatible Dataset
#'
#' Creates main processed dataset that's fully compatible with infill_soil_property()
#' Uses Module 8 utilities for data manipulation
#'
#' @param raw_data Original raw data
#' @param horizon_processing Horizon processing results
#' @param component_processing Component processing results
#' @param options Processing options
#' @param verbose Logical; provide progress messages
#'
#' @return Data frame compatible with existing infill workflow
#'
#' @export
create_infill_compatible_dataset <- function(raw_data,
                                             horizon_processing,
                                             component_processing,
                                             options,
                                             verbose = FALSE) {

  if (verbose) log_message("DEBUG", "Creating dataset compatible with infill_soil_property()", category = "DatasetCreation")

  # Start with raw data to preserve original structure
  processed_data <- raw_data

  # Apply unsuitable horizon detection using Module 8
  if (options$detect_unsuitable) {
    processed_data$unsuitable_horizon <- is_unsuitable(processed_data, hzname_col = "hzname")
  }

  # Apply data cleaning to preserve structure but improve quality
  if (options$advanced_cleaning) {
    property_columns <- identify_soil_property_columns_working(processed_data)

    for (prop_base in property_columns) {
      cleaning_result <- clean_property_data_ssurgo_compatible(
        processed_data, prop_base,
        generate_report = FALSE,
        verbose = FALSE
      )
      processed_data <- cleaning_result$data
    }
  }

  # Ensure essential columns exist and are properly typed
  processed_data <- ensure_essential_columns_working(processed_data)

  # Apply depth filtering if needed
  if ("hzdepb_r" %in% names(processed_data)) {
    max_depth <- options$max_depth %||% 250
    processed_data <- processed_data |>
      dplyr::filter(is.na(hzdepb_r) | hzdepb_r <= max_depth)
  }

  # Sort data for consistency with existing workflows
  if (all(c("cokey", "hzdept_r") %in% names(processed_data))) {
    processed_data <- processed_data |>
      dplyr::arrange(cokey, hzdept_r)
  }

  if (verbose) {
    log_message("DEBUG", paste("Final dataset:", nrow(processed_data), "rows"), category = "DatasetCreation")

    if ("unsuitable_horizon" %in% names(processed_data)) {
      unsuitable_count <- sum(processed_data$unsuitable_horizon, na.rm = TRUE)
      log_message("DEBUG", paste("Unsuitable horizons flagged:", unsuitable_count), category = "DatasetCreation")
    }
  }

  return(processed_data)
}

# ==============================================================================
# SSURGO-SPECIFIC CLEANING FUNCTIONS
# ==============================================================================

#' Clean Property Data (SSURGO Compatible)
#'
#' Enhanced version incorporating SSURGO-specific cleaning logic while using Module 8 utilities
#'
#' @param df Input data frame
#' @param property_name Name of the property to clean
#' @param validation_config Optional validation configuration
#' @param generate_report Whether to generate detailed quality report
#' @param verbose Logical; provide progress messages
#'
#' @return List containing cleaned data and optional quality report
#'
#' @seealso [advanced_string_parser_vectorized()], [vectorized_type_conversion()],
#'   and [apply_basic_range_limits()] - the string-parsing/range-limit helpers
#'   this function calls are defined in `data-infilling.R` (migrated from mod03),
#'   not in this file, after consolidating this file's near-duplicate
#'   `*_working()`/`*_ssurgo()` versions into those canonical implementations.
#' @export
clean_property_data_ssurgo_compatible <- function(df, property_name,
                                                  validation_config = NULL,
                                                  generate_report = TRUE,
                                                  verbose = FALSE) {

  start_time <- Sys.time()
  original_df <- df

  property_cols <- paste0(property_name, c("_r", "_l", "_h"))
  available_cols <- intersect(property_cols, names(df))

  if (length(available_cols) == 0) {
    if (verbose) log_message("DEBUG", paste("No property columns found for", property_name), category = "Cleaning")
    return(list(data = df))
  }

  # Initialize tracking
  cleaning_actions <- list()
  type_conversions <- list()
  outliers_detected <- list()
  parsing_results <- list()

  for (col in available_cols) {
    original_vals <- df[[col]]
    original_type <- class(original_vals)[1]

    # Advanced string parsing for character/factor data using Module 8 utilities
    if (is.character(original_vals) || is.factor(original_vals)) {
      parsing_result <- advanced_string_parser_vectorized(original_vals)
      df[[col]] <- parsing_result$values
      parsing_results[[col]] <- parsing_result$metadata
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

    # Statistical outlier detection using Module 8 (for _r columns)
    if (col == paste0(property_name, "_r")) {

      # Apply validation rules if provided
      if (!is.null(validation_config)) {
        property_ranges <- list()
        property_ranges[[property_name]] <- validation_config
        validation_result <- validate_numeric_ranges(df[col], property_ranges, action = "warn")
        if (length(validation_result$range_violations) > 0) {
          n_violations <- validation_result$range_violations[[property_name]]$n_violations
          cleaning_actions[[paste0(col, "_rule_violations")]] <- n_violations
        }
      }

      # Use Module 8 outlier detection
      outliers <- detect_outliers(df[[col]], method = "iqr", threshold = 3.0)
      df[[col]][outliers] <- NA
      outliers_detected[[col]] <- list(
        outliers = outliers,
        n_outliers = sum(outliers, na.rm = TRUE)
      )

      if (sum(outliers, na.rm = TRUE) > 0) {
        cleaning_actions[[paste0(col, "_statistical_outliers")]] <- sum(outliers, na.rm = TRUE)
      }

      # Apply basic range limits using SSURGO-specific ranges
      df[[col]] <- apply_basic_range_limits(df[[col]], property_name)
    }
  }

  # Generate quality report if requested
  if (generate_report) {
    end_time <- Sys.time()

    report <- list(
      property = property_name,
      timestamp = end_time,
      processing_time = difftime(end_time, start_time, units = "secs"),
      rows_processed = nrow(df),
      cleaning_actions = cleaning_actions,
      type_conversions = type_conversions,
      outliers_detected = outliers_detected,
      parsing_results = parsing_results
    )

    return(list(data = df, report = report))
  }

  return(list(data = df))
}

# ==============================================================================
# SSURGO-SPECIFIC UTILITY FUNCTIONS
# ==============================================================================
#
# NOTE: this file previously defined its own `advanced_string_parser_vectorized_working()`
# / `parse_single_string_advanced_working()` / `vectorized_type_conversion_working()` /
# `apply_basic_range_limits_ssurgo()`, near-duplicates of the canonical versions
# migrated from mod03 into data-infilling.R. They have been consolidated: this
# file now calls `advanced_string_parser_vectorized()`, `vectorized_type_conversion()`,
# and `apply_basic_range_limits()` (data-infilling.R) directly - see
# @seealso on clean_property_data_ssurgo_compatible() below.

#' Identify Soil Property Columns
#'
#' Finds base soil-property names (a known SSURGO property with an `_r`
#' representative-value column present in `df`).
#'
#' @param df Horizon data frame.
#' @return Character vector of property base names (e.g. `"sandtotal"`).
identify_soil_property_columns_working <- function(df) {
  # Identify base names of soil property columns (without _r, _l, _h suffixes)
  all_cols <- names(df)
  property_cols <- all_cols[grepl("_r$", all_cols)]
  property_bases <- gsub("_r$", "", property_cols)

  # Filter to known soil properties (SSURGO-specific)
  known_properties <- c("sandtotal", "claytotal", "silttotal", "dbovendry", "dbthirdbar",
                        "ph1to1h2o", "cec7", "om", "rfv", "wthirdbar", "wfifteenbar",
                        "awc", "ksat", "fragvol")

  return(intersect(property_bases, known_properties))
}

#' Calculate Property Completeness
#'
#' Fraction of non-`NA` `_r` values per known soil property in `df`.
#'
#' @param df Horizon data frame.
#' @return Named list of `property -> completeness fraction` (`[0, 1]`).
calculate_property_completeness_working <- function(df) {
  property_cols <- identify_soil_property_columns_working(df)
  completeness <- list()

  for (prop in property_cols) {
    r_col <- paste0(prop, "_r")
    if (r_col %in% names(df)) {
      completeness[[prop]] <- sum(!is.na(df[[r_col]])) / nrow(df)
    }
  }

  return(completeness)
}

#' Remove Invalid Horizon Records
#'
#' Drops horizon rows with missing/inconsistent depths (negative top depth,
#' bottom not greater than top, either depth beyond `max_depth`) or a
#' missing/empty `cokey`.
#'
#' @param df Horizon data frame.
#' @param max_depth Maximum plausible depth (cm).
#' @param verbose Logical; log the number of rows removed.
#' @return `df` with invalid rows removed.
remove_invalid_horizons_working_compatible <- function(df, max_depth = 250, verbose = FALSE) {
  initial_rows <- nrow(df)

  # Remove records with invalid depths
  if (all(c("hzdept_r", "hzdepb_r") %in% names(df))) {
    df <- df |>
      dplyr::filter(
        !is.na(hzdept_r),
        !is.na(hzdepb_r),
        hzdept_r >= 0,
        hzdepb_r > hzdept_r,
        hzdepb_r <= max_depth,
        hzdept_r < max_depth
      )
  }

  # Remove records with missing cokey
  if ("cokey" %in% names(df)) {
    df <- df |>
      dplyr::filter(!is.na(cokey), cokey != "")
  }

  if (verbose) {
    removed <- initial_rows - nrow(df)
    if (removed > 0) {
      log_message("DEBUG", paste("Removed", removed, "invalid horizon records"), category = "DataCleaning")
    }
  }

  return(df)
}

#' Remove Invalid Component Records
#'
#' Drops component rows with a missing/empty `cokey` or an out-of-range
#' `comppct_r` (component percentage, expected `[0, 100]`).
#'
#' @param df Component data frame.
#' @param verbose Logical; log the number of rows removed.
#' @return `df` with invalid rows removed.
remove_invalid_components_working <- function(df, verbose = FALSE) {
  initial_rows <- nrow(df)

  # Remove records with missing cokey
  if ("cokey" %in% names(df)) {
    df <- df |>
      dplyr::filter(!is.na(cokey), cokey != "")
  }

  # Remove records with invalid component percentages
  if ("comppct_r" %in% names(df)) {
    df <- df |>
      dplyr::filter(
        is.na(comppct_r) | (comppct_r >= 0 & comppct_r <= 100)
      )
  }

  if (verbose) {
    removed <- initial_rows - nrow(df)
    if (removed > 0) {
      log_message("DEBUG", paste("Removed", removed, "invalid component records"), category = "DataCleaning")
    }
  }

  return(df)
}

#' Calculate Basic Component Statistics
#'
#' Adds `is_major_component` (`comppct_r >= 15`) and, if `taxclname` is
#' present, `has_taxonomic_classification` flag columns.
#'
#' @param df Component data frame.
#' @param verbose Logical; log a completion message.
#' @return `df` with the derived columns added.
calculate_component_stats_working <- function(df, verbose = FALSE) {
  if (nrow(df) == 0) {
    return(df)
  }

  # Basic component statistics compatible with working code
  if ("comppct_r" %in% names(df)) {
    df$is_major_component <- ifelse(
      !is.na(df$comppct_r) & df$comppct_r >= 15,
      TRUE, FALSE
    )
  }

  # Add simple derived metrics
  if ("taxclname" %in% names(df)) {
    df$has_taxonomic_classification <- !is.na(df$taxclname)
  }

  if (verbose) {
    log_message("DEBUG", "Calculated basic component statistics", category = "ComponentProcessing")
  }

  return(df)
}

#' Calculate Derived Horizon Properties
#'
#' Fills in `hzthk_r` (thickness), `hz_midpoint`, and `awc_r` (available
#' water capacity, `wthirdbar_r - wfifteenbar_r`) when the source columns
#' they're derived from are present and the derived column is not already
#' present.
#'
#' @param df Horizon data frame.
#' @param verbose Logical; log a completion message.
#' @return `df` with the derived columns added.
calculate_derived_horizon_properties_working <- function(df, verbose = FALSE) {
  if (nrow(df) == 0) {
    return(df)
  }

  # Calculate horizon thickness if not present
  if (all(c("hzdept_r", "hzdepb_r") %in% names(df)) && !"hzthk_r" %in% names(df)) {
    df$hzthk_r <- df$hzdepb_r - df$hzdept_r
  }

  # Calculate horizon midpoint
  if (all(c("hzdept_r", "hzdepb_r") %in% names(df))) {
    df$hz_midpoint <- (df$hzdept_r + df$hzdepb_r) / 2
  }

  # Calculate available water capacity if water retention data available
  if (all(c("wthirdbar_r", "wfifteenbar_r") %in% names(df)) && !"awc_r" %in% names(df)) {
    df$awc_r <- df$wthirdbar_r - df$wfifteenbar_r
  }

  if (verbose) {
    log_message("DEBUG", "Calculated derived horizon properties", category = "HorizonProcessing")
  }

  return(df)
}

#' Ensure Essential Columns Are Present and Typed
#'
#' Guarantees `cokey` and `hzname` exist (synthesizing placeholder values if
#' missing) and coerces `hzdept_r`/`hzdepb_r` to numeric when present.
#'
#' @param df Horizon data frame.
#' @return `df` with the essential columns present/coerced.
ensure_essential_columns_working <- function(df) {
  # Ensure essential columns exist and have proper types

  # Ensure cokey exists
  if (!"cokey" %in% names(df)) {
    df$cokey <- paste0("missing_", seq_len(nrow(df)))
  }

  # Ensure hzname exists
  if (!"hzname" %in% names(df)) {
    df$hzname <- "Unknown"
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

#' Generate SSURGO Processing Quality Report
#'
#' Blends whichever quality signals are available (validation score, row
#' retention rate, horizon property completeness, component retention) into
#' a single weighted `[0, 1]` score and letter grade, renormalizing weights
#' over the signals actually present.
#'
#' @param original_data Data frame prior to processing.
#' @param processed_data Data frame after processing.
#' @param horizon_stats List possibly containing `property_completeness`.
#' @param component_stats List possibly containing `rows_retained`.
#' @param validation_results List possibly containing `overall_quality$score`.
#' @return List with `overall_quality_score`, `quality_grade`,
#'   `data_retention_rate`, a human-readable `processing_summary`, and the
#'   passed-through `horizon_processing`/`component_processing`/
#'   `validation_summary` inputs.
generate_processing_quality_report <- function(original_data, processed_data,
                                                horizon_stats, component_stats,
                                                validation_results) {
  data_retention_rate <- nrow(processed_data) / nrow(original_data)

  # Blend whichever quality signals are actually available into a single 0-1
  # score, rather than a hardcoded constant. Each component is itself already
  # a 0-1 fraction; weights are renormalized over the components present so
  # a missing signal (e.g. no property_completeness) doesn't just zero out
  # the score. NULL inputs are coalesced to NA_real_ explicitly (not left to
  # c()'s default behavior of silently dropping NULL elements, which would
  # misalign score_components against score_weights).
  na_if_null <- function(x) if (is.null(x)) NA_real_ else x

  score_components <- c(
    validation = na_if_null(validation_results$overall_quality$score),
    retention = min(data_retention_rate, 1),
    horizon_completeness = if (!is.null(horizon_stats$property_completeness) &&
                                length(horizon_stats$property_completeness) > 0) {
      mean(unlist(horizon_stats$property_completeness), na.rm = TRUE)
    } else {
      NA_real_
    },
    component_retention = na_if_null(component_stats$rows_retained)
  )
  score_weights <- c(validation = 0.4, retention = 0.25,
                      horizon_completeness = 0.25, component_retention = 0.10)

  available <- !is.na(score_components)
  overall_quality_score <- if (any(available)) {
    sum(score_components[available] * score_weights[available]) / sum(score_weights[available])
  } else {
    NA_real_
  }

  quality_grade <- if (is.na(overall_quality_score)) {
    "Unknown"
  } else if (overall_quality_score >= 0.9) {
    "Excellent"
  } else if (overall_quality_score >= 0.8) {
    "Good"
  } else if (overall_quality_score >= 0.7) {
    "Acceptable"
  } else if (overall_quality_score >= 0.6) {
    "Needs Improvement"
  } else {
    "Poor"
  }

  processing_summary <- sprintf(
    "Processed %d of %d rows (%.1f%% retained); overall quality score %s (%s)",
    nrow(processed_data), nrow(original_data), data_retention_rate * 100,
    if (is.na(overall_quality_score)) "NA" else sprintf("%.2f", overall_quality_score),
    quality_grade
  )

  return(list(
    overall_quality_score = overall_quality_score,
    quality_grade = quality_grade,
    data_retention_rate = data_retention_rate,
    processing_summary = processing_summary,
    horizon_processing = horizon_stats,
    component_processing = component_stats,
    validation_summary = validation_results,
    compatible_with_infill = TRUE
  ))
}

#' Calculate Horizon Quantiles and Probabilities by Mukey
#'
#' Calculates quantiles (5th, 50th, 95th percentiles) and prediction interval widths for soil
#' properties by mukey and depth, across a wide property basket (sand/silt/clay, bulk density,
#' water retention, RFV, pH, CEC, SOC, optional van Genuchten params), and, when texture data and
#' the optional `soiltexture` package are both available, the most probable USDA soil texture
#' class. Ported from `code/sim-functions.R:2728`.
#'
#' @param hz_data A data frame of simulated horizon data, with `mukey`, `hzdept_r`, `hzdepb_r`,
#'   and any of the recognized property columns (see the `prop_mapping`/`optional_props`
#'   translation table in the source).
#' @return A data frame with one row per `mukey`/depth, quantile statistics (`_05`/`_50`/`_95`/
#'   `_PIW90` suffixed columns), and - when texture data and `soiltexture` are available - the
#'   most probable texture class and its simulation-frequency probability.
#' @export
hz_quant_prob_mukey <- function(hz_data) {
  q <- c(0.05, 0.5, 0.95)

  # Start with base columns that we need
  base_cols <- c("mukey", "hzdept_r", "hzdepb_r")

  # Check which properties are available and build selection list
  prop_mapping <- list(
    "sand_total" = "sand",
    "silt_total" = "silt",
    "clay_total" = "clay",
    "db" = "Db",
    "wr_3b" = "wr3bar",
    "wr_15b" = "wr15bar",
    "rfv" = "rfv",
    "ph" = "ph",
    "cec" = "cec",
    "soc" = "soc"
  )

  # Optional properties that may not be present
  optional_props <- list(
    "theta_r" = "theta_r",
    "theta_s" = "theta_s",
    "alpha" = "alpha",
    "npar" = "npar",
    "ksat" = "ksat"
  )

  # Combine all possible properties
  all_props <- c(prop_mapping, optional_props)

  # Find which properties exist in the data
  available_props <- names(all_props)[names(all_props) %in% names(hz_data)]

  if (length(available_props) == 0) {
    stop("No valid soil properties found in the data")
  }

  # Create a renaming vector for the available properties (new_name = old_name)
  rename_vector <- character(0)
  for (col_name in available_props) {
    new_name <- all_props[[col_name]]
    rename_vector[new_name] <- col_name  # new_name = old_name format
  }

  # Select columns and rename them
  data_out <- hz_data |>
    dplyr::select(
      mukey = mukey,
      top = hzdept_r,
      bottom = hzdepb_r,
      dplyr::all_of(available_props)
    ) |>
    dplyr::rename(dplyr::all_of(rename_vector))

  # Get property names (excluding mukey, top, bottom)
  prop_names <- names(data_out)[!(names(data_out) %in% c("mukey", "top", "bottom"))]

  if (length(prop_names) == 0) {
    stop("No valid soil properties found in the data")
  }

  # PERF: previously three independent group_by(mukey, top)/summarize() passes over the same
  # data (one per quantile level) + three left_join()s to stitch them back together - replaced
  # with one grouped pass computing all three quantiles per group via across()'s multi-function
  # form, then a vectorized PIW90 (see PERFORMANCE_IMPROVEMENT_PLAN.md Tier 3). Column names
  # (`.names = "{.col}_{.fn}"`) match the "_05"/"_50"/"_95" suffixes the old code applied
  # afterward, so the final alphabetical column reorder below is unaffected.
  data_stats <- data_out |>
    dplyr::group_by(mukey, top) |>
    dplyr::summarize(
      dplyr::across(
        dplyr::all_of(prop_names),
        list(`05` = ~stats::quantile(.x, probs = q[1], na.rm = TRUE),
             `50` = ~stats::quantile(.x, probs = q[2], na.rm = TRUE),
             `95` = ~stats::quantile(.x, probs = q[3], na.rm = TRUE)),
        .names = "{.col}_{.fn}"
      ),
      .groups = "drop"
    ) |>
    as.data.frame()

  # Calculate prediction interval width (PIW90) for each property, from the columns just computed.
  for (prop in prop_names) {
    data_stats[[paste0(prop, "_PIW90")]] <-
      data_stats[[paste0(prop, "_95")]] - data_stats[[paste0(prop, "_05")]]
  }

  # Reorder columns alphabetically (excluding mukey and top)
  other_cols <- names(data_stats)[!(names(data_stats) %in% c("mukey", "top"))]
  data_stats <- data_stats |>
    dplyr::select(mukey, top, dplyr::all_of(other_cols[order(other_cols)]))

  # Add bottom depth information
  depth_info <- hz_data |>
    dplyr::select(mukey, top = hzdept_r, bottom = hzdepb_r) |>
    dplyr::distinct()
  data_stats <- depth_info |>
    dplyr::left_join(data_stats, by = c("mukey", "top"))

  # Calculate texture class probabilities if texture data is available
  if (all(c("sand", "silt", "clay") %in% names(data_out))) {

    txt_out <- data_out |>
      dplyr::select(mukey, top, SAND = sand, SILT = silt, CLAY = clay)

    # Check if soiltexture package is available
    if (requireNamespace("soiltexture", quietly = TRUE)) {
      txt_out$txt_class <- soiltexture::TT.points.in.classes(
        tri.data = txt_out,
        class.sys = "USDA-NCSS.TT",
        PiC.type = "t"
      )

      # Calculate texture class probabilities
      n_sims <- length(unique(hz_data$simulation_number))
      if (is.na(n_sims) || n_sims == 0) n_sims <- 10  # fallback

      txt_class_prob <- txt_out |>
        dplyr::group_by(mukey, top) |>
        dplyr::count(txt_class) |>
        dplyr::mutate(prob = n / n_sims) |>
        dplyr::ungroup() |>
        dplyr::group_by(mukey, top) |>
        dplyr::slice_max(prob, n = 1) |>
        dplyr::ungroup() |>
        as.data.frame()

      # Add texture class information
      txt_class_final <- depth_info |>
        dplyr::left_join(txt_class_prob |> dplyr::select(-n), by = c("mukey", "top"))

      data_stats <- data_stats |>
        dplyr::left_join(
          txt_class_final |> dplyr::select(-bottom),
          by = c("mukey", "top")
        )
    } else {
      warning("soiltexture package not available. Texture class probabilities will not be calculated.")
    }
  }

  return(data_stats)
}
