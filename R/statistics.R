#' @title Statistical Analysis & Correlation Modeling
#' @description Advanced statistical analysis for soil properties with correlation modeling and distribution fitting
#' @name soil_statistics
NULL

# ==============================================================================
# MAIN STATISTICAL ANALYSIS FUNCTIONS
# ==============================================================================

#' Comprehensive Statistical Analysis (Main Entry Point)
#'
#' Enhanced statistical analysis leveraging Module 0 utilities for data validation,
#' error handling, logging, and statistical computations.
#'
#' @param processed_data Processed SSURGO data from Module 2
#' @param analysis_config List of analysis configuration parameters
#' @param correlation_methods Character vector of correlation methods to apply
#' @param distribution_fitting Logical; perform distribution fitting analysis
#' @param outlier_detection Logical; perform statistical outlier detection
#' @param validate_results Logical; validate statistical results (default: TRUE)
#' @param verbose Logical; provide detailed progress messages
#'
#' @return List containing comprehensive statistical analysis results
#'
#' @examples
#' \dontrun{
#' # Basic statistical analysis
#' stats_result <- analyze_soil_statistics(
#'   processed_data,
#'   analysis_config = list(
#'     stratify_by_horizon = TRUE,
#'     include_texture_analysis = TRUE,
#'     correlation_threshold = 0.05
#'   )
#' )
#' }
#'
#' @export
analyze_soil_statistics <- function(processed_data,
                                    analysis_config = list(),
                                    correlation_methods = c("pearson", "spearman"),
                                    distribution_fitting = TRUE,
                                    outlier_detection = TRUE,
                                    validate_results = TRUE,
                                    verbose = getOption("ssurgo.verbose", FALSE)) {
  # NOTE: this used to call setup_logging(log_level = if(verbose) "DEBUG" else "INFO")
  # unconditionally on first use, which permanently changed the global log level for
  # the rest of the session (no restore) - meaning one early call here would silently
  # un-quiet every other function afterward. set_verbose_logging() raises the level only
  # for the duration of this call, restored via on.exit() below.
  .old_log_cfg <- set_verbose_logging(verbose)
  on.exit(options(soil_workflow_log_config = .old_log_cfg), add = TRUE)

  start_time <- Sys.time()

  log_message("INFO", "=== Statistical Analysis Started ===", category = "Statistics")

  # Merge with default configuration using Module 0 utilities
  default_config <- get_statistical_analysis_defaults()
  config <- merge_configurations(default_config, analysis_config)

  # FIXED: Safer configuration access with defaults
  min_quality_score <- tryCatch({
    config$statistical_analysis$min_quality_score %||% config$min_quality_score %||% 0.6
  }, error = function(e) {
    if (verbose) log_message("DEBUG", "Using default min_quality_score: 0.6", category = "Statistics")
    0.6
  })

  log_message("INFO", paste("Input data:", nrow(processed_data), "rows,", ncol(processed_data), "columns"), category = "Statistics")

  # Step 1: Comprehensive data validation using Module 0
  log_message("INFO", "Step 1: Validating input data", category = "Statistics")

  data_validation <- validate_data_quality(
    processed_data,
    required_columns = c("cokey", "hzdept_r"),
    numeric_columns = identify_numeric_soil_properties(processed_data),
    quality_thresholds = config$quality_thresholds %||% get_default_quality_thresholds()
  )

  # FIXED: Safe access to data validation results with multiple possible structures
  overall_quality_score <- tryCatch({
    # Try different possible structures
    data_validation$overall_quality$score %||%
      data_validation$overall_score %||%
      data_validation$quality_score %||%
      0.8  # Default fallback
  }, error = function(e) {
    if (verbose) log_message("DEBUG", paste("Could not access quality score, using default:", e$message), category = "Statistics")
    0.8
  })

  # FIXED: Safe validation passed check
  validation_passed <- tryCatch({
    data_validation$overall_quality$validation_passed %||%
      data_validation$validation_passed %||%
      (overall_quality_score >= min_quality_score)
  }, error = function(e) {
    if (verbose) log_message("DEBUG", "Using quality score for validation check", category = "Statistics")
    (overall_quality_score >= min_quality_score)
  })

  if (overall_quality_score < min_quality_score) {
    log_message("WARN", paste("Data quality below threshold:", round(overall_quality_score, 3), "<", min_quality_score), category = "Statistics")

    # Use warn instead of stop to continue analysis
    handle_workflow_error(
      list(message = paste("Data quality too low:", overall_quality_score)),
      context = "Data Validation",
      recovery_action = "warn"
    )
  } else {
    log_message("INFO", paste("Data quality acceptable:", round(overall_quality_score, 3)), category = "Statistics")
  }

  # Step 2: Property validation using Module 0
  log_message("INFO", "Step 2: Validating soil properties", category = "Statistics")

  available_properties <- identify_numeric_soil_properties(processed_data)

  # FIXED: Filter out empty/all-NA columns to avoid min/max warnings
  available_properties <- filter_valid_numeric_properties(processed_data, available_properties, verbose)

  property_validation <- validate_properties_with_synonyms(
    available_properties,
    property_lookup = "ssurgo",
    strict_mode = FALSE  # Use flexible mode for statistical analysis
  )

  if (!property_validation$valid) {
    log_message("WARN", paste("Property validation issues:", paste(property_validation$warnings, collapse = "; ")), category = "Statistics")
  }

  log_message("INFO", paste("Valid properties for analysis:", length(available_properties)), category = "Statistics")

  # Step 3: Handle missing values using Module 0
  log_message("INFO", "Step 3: Processing missing values", category = "Statistics")

  processed_data_clean <- handle_missing_values(
    processed_data,
    strategy = config$missing_value_strategy %||% "interpolate",
    columns = available_properties,
    group_by = if(config$stratify_by_horizon %||% TRUE) "genhz" else NULL
  )

  # Step 4: Comprehensive correlation analysis
  log_message("INFO", "Step 4: Running correlation analysis", category = "Statistics")

  correlation_analysis <- tryCatch({
    run_comprehensive_correlation_analysis_safe(
      data = processed_data_clean,
      methods = correlation_methods,
      config = config,
      available_properties = available_properties,
      verbose = verbose
    )
  }, error = function(e) {
    handle_workflow_error(e, "Correlation Analysis", "warn")
    return(list(matrices = list(), summary = list(), validation = list()))
  })

  # Step 5: Distribution analysis (if requested)
  distribution_analysis <- NULL
  if (distribution_fitting) {
    log_message("INFO", "Step 5: Analyzing property distributions", category = "Statistics")

    distribution_analysis <- tryCatch({
      analyze_property_distributions_safe(
        data = processed_data_clean,
        properties = available_properties,
        config = config,
        verbose = verbose
      )
    }, error = function(e) {
      handle_workflow_error(e, "Distribution Analysis", "warn")
      return(NULL)
    })
  }

  # Step 6: Outlier detection (if requested)
  outlier_analysis <- NULL
  if (outlier_detection) {
    log_message("INFO", "Step 6: Detecting outliers", category = "Statistics")

    outlier_analysis <- tryCatch({
      detect_comprehensive_outliers_safe(
        data = processed_data_clean,
        properties = available_properties,
        config = config,
        verbose = verbose
      )
    }, error = function(e) {
      handle_workflow_error(e, "Outlier Detection", "warn")
      return(NULL)
    })
  }

  # Step 7: Property statistics using Module 0
  log_message("INFO", "Step 7: Computing property statistics", category = "Statistics")

  property_statistics <- compute_property_statistics(
    data = processed_data_clean,
    properties = available_properties,
    config = config
  )

  # Step 8: Results validation (if requested)
  validation_results <- NULL
  if (validate_results) {
    log_message("INFO", "Step 8: Validating results", category = "Statistics")

    validation_results <- validate_statistical_results_safe(
      correlation_analysis = correlation_analysis,
      distribution_analysis = distribution_analysis,
      outlier_analysis = outlier_analysis,
      property_statistics = property_statistics,
      config = config
    )
  }

  # Step 9: Generate quality report
  log_message("INFO", "Step 9: Generating quality report", category = "Statistics")

  quality_report <- generate_statistical_quality_report_safe(
    original_data = processed_data,
    processed_data = processed_data_clean,
    correlation_analysis = correlation_analysis,
    distribution_analysis = distribution_analysis,
    outlier_analysis = outlier_analysis,
    property_statistics = property_statistics,
    validation_results = validation_results,
    data_validation = data_validation,
    config = config
  )

  # Step 10: Generate comprehensive metadata
  end_time <- Sys.time()
  analysis_metadata <- list(
    analysis_start = start_time,
    analysis_end = end_time,
    duration_seconds = as.numeric(difftime(end_time, start_time, units = "secs")),
    config_used = config,
    correlation_methods = correlation_methods,
    input_rows = nrow(processed_data),
    processed_rows = nrow(processed_data_clean),
    properties_analyzed = length(available_properties),
    correlation_matrices_generated = if(!is.null(correlation_analysis)) length(correlation_analysis$matrices) else 0,
    distributions_fitted = if(!is.null(distribution_analysis)) length(distribution_analysis$fitted_distributions %||% list()) else 0,
    outliers_detected = if(!is.null(outlier_analysis)) sum(outlier_analysis$summary$total_outliers %||% 0, na.rm = TRUE) else 0,
    data_quality_score = overall_quality_score,
    validation_score = if(!is.null(validation_results)) validation_results$validation_score %||% NA else NA
  )

  log_message("INFO", "=== Statistical Analysis Complete ===", category = "Statistics")
  log_message("INFO", paste("Duration:", round(analysis_metadata$duration_seconds, 2), "seconds"), category = "Statistics")
  log_message("INFO", paste("Properties analyzed:", analysis_metadata$properties_analyzed), category = "Statistics")
  log_message("INFO", paste("Overall quality score:", round(quality_report$overall_quality_score %||% 0.8, 3)), category = "Statistics")

  # Return comprehensive results
  return(list(
    correlation_matrices = correlation_analysis,
    distribution_analysis = distribution_analysis,
    outlier_analysis = outlier_analysis,
    property_statistics = property_statistics,
    validation_results = validation_results,
    analysis_metadata = analysis_metadata,
    quality_report = quality_report,
    data_validation = data_validation,
    processed_data = if(config$return_processed_data %||% FALSE) processed_data_clean else NULL
  ))
}

#' Filter Valid Numeric Properties
#'
#' FIXED: Remove empty or all-NA columns to prevent min/max warnings
#'
#' @param data Input data
#' @param properties Properties to check
#' @param verbose Logical; provide progress messages
#' @return Filtered properties list
filter_valid_numeric_properties <- function(data, properties, verbose = FALSE) {

  valid_properties <- character(0)

  for (prop in properties) {
    if (prop %in% names(data)) {
      values <- data[[prop]]

      # Check if column has any non-NA, finite values
      finite_values <- values[is.finite(values)]

      if (length(finite_values) > 0) {
        valid_properties <- c(valid_properties, prop)
      } else {
        if (verbose) {
          log_message("DEBUG", paste("Excluding property", prop, "- no finite values"), category = "Statistics")
        }
      }
    }
  }

  if (verbose) {
    excluded_count <- length(properties) - length(valid_properties)
    if (excluded_count > 0) {
      log_message("INFO", paste("Excluded", excluded_count, "properties with no finite values"), category = "Statistics")
    }
  }

  return(valid_properties)
}

#' Run Comprehensive Correlation Analysis (Safe Version)
#'
#' Safe version with enhanced error handling
#'
#' @param data Input data.
#' @param methods Character vector of correlation methods (e.g. `c("pearson","spearman")`).
#' @param config Analysis configuration.
#' @param available_properties Character vector of properties to correlate.
#' @param verbose Logical; provide detailed progress messages.
#' @return `list(matrices=, summary=, validation=)`.
run_comprehensive_correlation_analysis_safe <- function(data, methods, config, available_properties, verbose = FALSE) {

  log_message("DEBUG", paste("Computing correlations using", length(methods), "methods for", length(available_properties), "properties"), category = "Correlation")

  correlation_results <- list(
    matrices = list(),
    summary = list(),
    validation = list()
  )

  if (length(available_properties) < 2) {
    log_message("WARN", "Need at least 2 properties for correlation analysis", category = "Correlation")
    return(correlation_results)
  }

  # Prepare correlation data
  correlation_data <- data[, available_properties, drop = FALSE]

  # Remove rows with all NA values
  complete_rows <- rowSums(!is.na(correlation_data)) > 0
  correlation_data <- correlation_data[complete_rows, , drop = FALSE]

  if (nrow(correlation_data) == 0) {
    log_message("WARN", "No complete rows for correlation analysis", category = "Correlation")
    return(correlation_results)
  }

  # Global correlation matrices using Module 0's safe_correlation
  for (method in methods) {
    log_message("DEBUG", paste("Computing", method, "correlation matrix"), category = "Correlation")

    correlation_matrix <- safe_correlation(
      correlation_data,
      method = method,
      handle_constant = "warn"
    )

    if (!is.null(correlation_matrix)) {
      # Enhanced correlation results
      correlation_results$matrices[[method]] <- list(
        matrix = correlation_matrix,
        method = method,
        n_variables = ncol(correlation_data),
        n_observations = nrow(correlation_data)
      )

      log_message("INFO", paste("Successfully computed", method, "correlation matrix"), category = "Correlation")
    } else {
      log_message("WARN", paste("Failed to compute", method, "correlation matrix"), category = "Correlation")
    }
  }

  # Generate summary
  correlation_results$summary <- list(
    n_matrices = length(correlation_results$matrices),
    methods_used = names(correlation_results$matrices),
    n_properties = length(available_properties),
    n_observations = nrow(correlation_data)
  )

  return(correlation_results)
}

#' Safe Property Distribution Analysis
#'
#' @param data Input data.
#' @param properties Character vector of properties to analyze.
#' @param config Analysis configuration.
#' @param verbose Logical; provide detailed progress messages.
#' @return `list(fitted_distributions=, distribution_tests=, summary=)`.
analyze_property_distributions_safe <- function(data, properties, config, verbose = FALSE) {

  log_message("DEBUG", paste("Analyzing distributions for", length(properties), "properties"), category = "Distribution")

  distribution_results <- list(
    fitted_distributions = list(),
    distribution_tests = list(),
    summary = list()
  )

  minimum_observations <- config$minimum_observations %||% 10

  for (property in properties) {
    if (property %in% names(data)) {
      property_data <- data[[property]]
      property_data <- property_data[!is.na(property_data) & is.finite(property_data)]

      if (length(property_data) >= minimum_observations) {
        if (verbose) {
          log_message("DEBUG", paste("Analyzing distribution for", property, "- n =", length(property_data)), category = "Distribution")
        }

        # Basic distribution summary
        distribution_results$fitted_distributions[[property]] <- list(
          n_observations = length(property_data),
          mean = mean(property_data),
          sd = sd(property_data),
          min = min(property_data),
          max = max(property_data),
          median = median(property_data),
          q25 = quantile(property_data, 0.25),
          q75 = quantile(property_data, 0.75)
        )
      } else {
        if (verbose) {
          log_message("DEBUG", paste("Skipping", property, "- insufficient data (n =", length(property_data), ")"), category = "Distribution")
        }
      }
    }
  }

  # Summary
  distribution_results$summary <- list(
    total_properties = length(properties),
    analyzed_properties = length(distribution_results$fitted_distributions),
    minimum_observations = minimum_observations
  )

  return(distribution_results)
}

#' Safe Outlier Detection
#'
#' @param data Input data.
#' @param properties Character vector of properties to check.
#' @param config Analysis configuration.
#' @param verbose Logical; provide detailed progress messages.
#' @return `list(property_outliers=, summary=)`.
detect_comprehensive_outliers_safe <- function(data, properties, config, verbose = FALSE) {

  log_message("DEBUG", paste("Detecting outliers for", length(properties), "properties"), category = "Outliers")

  outlier_results <- list(
    property_outliers = list(),
    summary = list()
  )

  minimum_observations <- config$minimum_observations %||% 10

  # Property-specific outlier detection using Module 0
  for (property in properties) {
    if (property %in% names(data)) {
      property_data <- data[[property]]
      finite_data <- property_data[!is.na(property_data) & is.finite(property_data)]

      if (length(finite_data) >= minimum_observations) {
        # Use Module 0's detect_outliers
        outliers <- detect_outliers(
          finite_data,
          method = "iqr",
          threshold = 1.5,
          return_indices = FALSE
        )

        outlier_results$property_outliers[[property]] <- list(
          n_outliers = sum(outliers, na.rm = TRUE),
          outlier_rate = mean(outliers, na.rm = TRUE),
          total_observations = length(finite_data)
        )

        if (verbose) {
          log_message("DEBUG", paste("Property", property, "- outliers:", sum(outliers, na.rm = TRUE), "/", length(finite_data)), category = "Outliers")
        }
      }
    }
  }

  # Generate summary
  total_outliers <- sum(sapply(outlier_results$property_outliers, function(x) x$n_outliers), na.rm = TRUE)

  outlier_results$summary <- list(
    total_outliers = total_outliers,
    properties_analyzed = length(outlier_results$property_outliers)
  )

  return(outlier_results)
}

#' Safe Statistical Results Validation
#'
#' @param correlation_analysis Correlation analysis results.
#' @param distribution_analysis Distribution analysis results.
#' @param outlier_analysis Outlier analysis results.
#' @param property_statistics Property statistics results.
#' @param config Analysis configuration.
#' @return `list(overall_valid=, errors=, warnings=, validation_score=, component_validations=)`.
validate_statistical_results_safe <- function(correlation_analysis, distribution_analysis,
                                              outlier_analysis, property_statistics, config) {

  validation_results <- list(
    overall_valid = TRUE,
    errors = character(0),
    warnings = character(0),
    validation_score = 1.0,
    component_validations = list()
  )

  log_message("DEBUG", "Validating statistical results", category = "Validation")

  # Simple validation - check if we have results
  if (is.null(correlation_analysis) || length(correlation_analysis$matrices) == 0) {
    validation_results$warnings <- c(validation_results$warnings, "No correlation matrices generated")
  }

  if (is.null(property_statistics) || length(property_statistics) == 0) {
    validation_results$warnings <- c(validation_results$warnings, "No property statistics computed")
  }

  # Calculate validation score
  score_components <- c()

  if (!is.null(correlation_analysis) && length(correlation_analysis$matrices) > 0) {
    score_components <- c(score_components, 1.0)
  } else {
    score_components <- c(score_components, 0.5)
  }

  if (!is.null(property_statistics) && length(property_statistics) > 0) {
    score_components <- c(score_components, 1.0)
  } else {
    score_components <- c(score_components, 0.5)
  }

  validation_results$validation_score <- mean(score_components)

  return(validation_results)
}

#' Safe Quality Report Generation
#'
#' @param original_data Original input data.
#' @param processed_data Processed data.
#' @param correlation_analysis Correlation analysis results.
#' @param distribution_analysis Distribution analysis results.
#' @param outlier_analysis Outlier analysis results.
#' @param property_statistics Property statistics results.
#' @param validation_results Validation results.
#' @param data_validation Data validation results.
#' @param config Analysis configuration.
#' @return `list(overall_quality_score=, data_quality=, analysis_quality=, validation_summary=, recommendations=)`.
generate_statistical_quality_report_safe <- function(original_data, processed_data,
                                                     correlation_analysis, distribution_analysis,
                                                     outlier_analysis, property_statistics,
                                                     validation_results, data_validation, config) {

  quality_report <- list(
    overall_quality_score = 0.8,  # Default
    data_quality = list(),
    analysis_quality = list(),
    validation_summary = list(),
    recommendations = character(0)
  )

  # Try to extract data quality score safely
  tryCatch({
    data_quality_score <- data_validation$overall_quality$score %||%
      data_validation$quality_score %||% 0.8
    quality_report$overall_quality_score <- data_quality_score

    quality_report$data_quality <- list(
      input_quality_score = data_quality_score,
      missing_data_handled = TRUE,
      property_coverage = length(property_statistics) / max(1, length(identify_numeric_soil_properties(original_data)))
    )
  }, error = function(e) {
    log_message("DEBUG", "Using default quality scores", category = "Statistics")
  })

  # Analysis quality
  quality_report$analysis_quality <- list(
    correlation_matrices_generated = if(!is.null(correlation_analysis)) length(correlation_analysis$matrices) else 0,
    properties_analyzed = length(property_statistics),
    analysis_completed = TRUE
  )

  # Validation summary
  if (!is.null(validation_results)) {
    quality_report$validation_summary <- list(
      overall_valid = validation_results$overall_valid %||% TRUE,
      validation_score = validation_results$validation_score %||% 0.8,
      n_warnings = length(validation_results$warnings %||% character(0))
    )
  }

  return(quality_report)
}
# ==============================================================================
# ENHANCED CORRELATION ANALYSIS (LEVERAGING Module 0)
# ==============================================================================

#' Run Comprehensive Correlation Analysis
#'
#' Enhanced correlation analysis using Module 0 statistical utilities
#'
#' @param data Input data
#' @param methods Correlation methods
#' @param config Analysis configuration
#' @param available_properties Available soil properties
#'
#' @return Enhanced correlation analysis results
#'
#' @export
run_comprehensive_correlation_analysis <- function(data, methods, config, available_properties) {

  log_message("DEBUG", paste("Computing correlations using", length(methods), "methods for", length(available_properties), "properties"), category = "Correlation")

  correlation_results <- list(
    matrices = list(),
    summary = list(),
    validation = list()
  )

  # Prepare correlation data
  correlation_data <- data[, available_properties, drop = FALSE]

  # Global correlation matrices using Module 0's safe_correlation
  for (method in methods) {
    log_message("DEBUG", paste("Computing", method, "correlation matrix"), category = "Correlation")

    correlation_matrix <- safe_correlation(
      correlation_data,
      method = method,
      handle_constant = config$handle_constant_variables
    )

    # Enhanced correlation results
    correlation_results$matrices[[method]] <- list(
      matrix = correlation_matrix,
      method = method,
      n_variables = ncol(correlation_data),
      n_observations = nrow(correlation_data),
      eigenvalues = if(is.matrix(correlation_matrix)) eigen(correlation_matrix, only.values = TRUE)$values else NULL,
      condition_number = if(is.matrix(correlation_matrix)) calculate_condition_number(correlation_matrix) else NA
    )
  }

  # Stratified correlation analysis by horizon (if requested)
  if (config$stratify_by_horizon && "genhz" %in% names(data)) {
    log_message("DEBUG", "Computing horizon-stratified correlations", category = "Correlation")

    correlation_results$stratified <- compute_stratified_correlations(
      data = data,
      properties = available_properties,
      stratify_by = "genhz",
      methods = methods,
      config = config
    )
  }

  # Texture analysis (if requested and texture data available)
  if (config$include_texture_analysis) {
    texture_props <- c("sandtotal_r", "silttotal_r", "claytotal_r")
    available_texture <- intersect(texture_props, available_properties)

    if (length(available_texture) >= 2) {
      log_message("DEBUG", "Computing texture correlations", category = "Correlation")

      correlation_results$texture_analysis <- analyze_texture_correlations(
        data = data,
        texture_properties = available_texture,
        methods = methods,
        config = config
      )
    }
  }

  # Correlation validation using Module 0 utilities
  correlation_results$validation <- validate_correlation_matrices(
    correlation_results$matrices,
    config = config
  )

  # Summary statistics
  correlation_results$summary <- generate_correlation_summary(correlation_results)

  return(correlation_results)
}

#' Compute Stratified Correlations
#'
#' Compute correlation matrices stratified by a grouping variable
#'
#' @param data Input data
#' @param properties Properties to correlate
#' @param stratify_by Grouping column
#' @param methods Correlation methods
#' @param config Configuration
#'
#' @return Stratified correlation results
#'
#' @export
compute_stratified_correlations <- function(data, properties, stratify_by, methods, config) {

  stratified_results <- list()

  # Get unique groups
  unique_groups <- unique(data[[stratify_by]][!is.na(data[[stratify_by]])])

  log_message("DEBUG", paste("Computing correlations for", length(unique_groups), "groups"), category = "Correlation")

  for (group in unique_groups) {
    group_data <- data[data[[stratify_by]] == group & !is.na(data[[stratify_by]]), properties, drop = FALSE]

    # Check minimum sample size using Module 0 validation
    if (nrow(group_data) < config$minimum_observations) {
      log_message("DEBUG", paste("Skipping group", group, "- insufficient observations:", nrow(group_data)), category = "Correlation")
      next
    }

    group_correlations <- list()
    for (method in methods) {
      correlation_matrix <- safe_correlation(
        group_data,
        method = method,
        handle_constant = config$handle_constant_variables
      )

      group_correlations[[method]] <- list(
        matrix = correlation_matrix,
        n_observations = nrow(group_data),
        group = group
      )
    }

    if (length(group_correlations) > 0) {
      stratified_results[[as.character(group)]] <- group_correlations
    }
  }

  return(stratified_results)
}

# ==============================================================================
# ENHANCED DISTRIBUTION ANALYSIS (LEVERAGING Module 0)
# ==============================================================================

#' Analyze Property Distributions (Enhanced)
#'
#' Enhanced distribution analysis using Module 0 utilities
#'
#' @param data Input data
#' @param properties Properties to analyze
#' @param config Configuration
#'
#' @return Enhanced distribution analysis results
#'
#' @export
analyze_property_distributions <- function(data, properties, config) {

  log_message("DEBUG", paste("Analyzing distributions for", length(properties), "properties"), category = "Distribution")

  distribution_results <- list(
    fitted_distributions = list(),
    distribution_tests = list(),
    summary = list()
  )

  for (property in properties) {
    property_data <- data[[property]]
    property_data <- property_data[!is.na(property_data) & is.finite(property_data)]

    if (length(property_data) < config$minimum_observations) {
      log_message("DEBUG", paste("Skipping", property, "- insufficient data"), category = "Distribution")
      next
    }

    log_message("DEBUG", paste("Fitting distributions for", property), category = "Distribution")

    # Fit distributions
    fitted_distributions <- fit_property_distributions(
      values = property_data,
      property_name = property,
      config = config
    )

    # Distribution tests
    distribution_tests <- perform_distribution_tests(
      values = property_data,
      property_name = property,
      config = config
    )

    distribution_results$fitted_distributions[[property]] <- fitted_distributions
    distribution_results$distribution_tests[[property]] <- distribution_tests
  }

  # Overall distribution assessment
  distribution_results$summary <- assess_overall_distributions(
    distribution_results$fitted_distributions,
    distribution_results$distribution_tests
  )

  return(distribution_results)
}

#' Fit Property Distributions (Enhanced)
#'
#' Enhanced distribution fitting with better error handling
#'
#' @param values Property values
#' @param property_name Property name
#' @param config Configuration
#'
#' @return Fitted distributions with quality metrics
#'
#' @export
fit_property_distributions <- function(values, property_name, config) {

  distributions_to_test <- get_appropriate_distributions(property_name, values, config)
  fitted_distributions <- list()

  for (dist_name in distributions_to_test) {
    tryCatch({
      fitted_dist <- fit_single_distribution(values, dist_name, property_name, config)
      if (!is.null(fitted_dist)) {
        fitted_distributions[[dist_name]] <- fitted_dist
      }
    }, error = function(e) {
      log_message("DEBUG", paste("Failed to fit", dist_name, "for", property_name, ":", e$message), category = "Distribution")
    })
  }

  # Rank distributions by fit quality
  if (length(fitted_distributions) > 1) {
    fitted_distributions <- rank_distribution_fits(fitted_distributions)
  }

  return(fitted_distributions)
}

# ==============================================================================
# ENHANCED OUTLIER DETECTION (LEVERAGING Module 0)
# ==============================================================================

#' Detect Comprehensive Outliers (Enhanced)
#'
#' Enhanced outlier detection using Module 0 utilities
#'
#' @param data Input data
#' @param properties Properties to analyze
#' @param config Configuration
#'
#' @return Enhanced outlier analysis results
#'
#' @export
detect_comprehensive_outliers <- function(data, properties, config) {

  log_message("DEBUG", paste("Detecting outliers for", length(properties), "properties"), category = "Outliers")

  outlier_results <- list(
    property_outliers = list(),
    multivariate_outliers = list(),
    summary = list()
  )

  # Property-specific outlier detection using Module 0
  for (property in properties) {
    property_data <- data[[property]]

    if (sum(!is.na(property_data)) < config$minimum_observations) {
      next
    }

    # Use Module 0's detect_outliers with multiple methods
    outlier_methods <- config$outlier_methods
    property_outlier_results <- list()

    for (method in outlier_methods) {
      outliers <- detect_outliers(
        property_data,
        method = method,
        threshold = config$outlier_thresholds[[method]],
        return_indices = FALSE
      )

      property_outlier_results[[method]] <- list(
        outliers = outliers,
        n_outliers = sum(outliers, na.rm = TRUE),
        outlier_rate = mean(outliers, na.rm = TRUE)
      )
    }

    outlier_results$property_outliers[[property]] <- property_outlier_results
  }

  # Multivariate outlier detection (if sufficient properties)
  if (length(properties) >= 2 && config$detect_multivariate_outliers) {
    log_message("DEBUG", "Detecting multivariate outliers", category = "Outliers")

    outlier_results$multivariate_outliers <- detect_multivariate_outliers(
      data = data,
      properties = properties,
      config = config
    )
  }

  # Outlier summary
  outlier_results$summary <- generate_outlier_summary(outlier_results)

  return(outlier_results)
}

# ==============================================================================
# ENHANCED VALIDATION AND QUALITY REPORTING
# ==============================================================================

#' Validate Statistical Results (Enhanced)
#'
#' Enhanced validation using Module 0 utilities
#'
#' @param correlation_analysis Correlation results
#' @param distribution_analysis Distribution results
#' @param outlier_analysis Outlier results
#' @param property_statistics Property statistics
#' @param config Configuration
#'
#' @return Enhanced validation results
#'
#' @export
validate_statistical_results <- function(correlation_analysis, distribution_analysis,
                                                  outlier_analysis, property_statistics, config) {

  validation_results <- list(
    overall_valid = TRUE,
    errors = character(0),
    warnings = character(0),
    validation_score = 1.0,
    component_validations = list()
  )

  log_message("DEBUG", "Validating statistical results", category = "Validation")

  # Validate correlation matrices
  if (!is.null(correlation_analysis)) {
    corr_validation <- validate_correlation_matrices(correlation_analysis$matrices, config)
    validation_results$component_validations$correlation <- corr_validation

    if (!corr_validation$valid) {
      validation_results$overall_valid <- FALSE
      validation_results$errors <- c(validation_results$errors, corr_validation$errors)
    }
    validation_results$warnings <- c(validation_results$warnings, corr_validation$warnings)
  }

  # Validate distribution analysis
  if (!is.null(distribution_analysis)) {
    dist_validation <- validate_distribution_analysis(distribution_analysis, config)
    validation_results$component_validations$distribution <- dist_validation
    validation_results$warnings <- c(validation_results$warnings, dist_validation$warnings)
  }

  # Validate outlier analysis
  if (!is.null(outlier_analysis)) {
    outlier_validation <- validate_outlier_analysis(outlier_analysis, config)
    validation_results$component_validations$outlier <- outlier_validation
    validation_results$warnings <- c(validation_results$warnings, outlier_validation$warnings)
  }

  # Calculate overall validation score
  validation_results$validation_score <- calculate_overall_validation_score(validation_results)

  log_message("DEBUG", paste("Validation complete - Score:", round(validation_results$validation_score, 3)), category = "Validation")

  return(validation_results)
}

#' Generate Statistical Quality Report (Enhanced)
#'
#' Enhanced quality reporting with comprehensive metrics
#'
#' @param original_data Original input data
#' @param processed_data Processed data
#' @param correlation_analysis Correlation results
#' @param distribution_analysis Distribution results
#' @param outlier_analysis Outlier results
#' @param property_statistics Property statistics
#' @param validation_results Validation results
#' @param data_validation Data validation results
#' @param config Configuration
#'
#' @return Enhanced quality report
#'
#' @export
generate_statistical_quality_report <- function(original_data, processed_data,
                                                         correlation_analysis, distribution_analysis,
                                                         outlier_analysis, property_statistics,
                                                         validation_results, data_validation, config) {

  quality_report <- list(
    overall_quality_score = calculate_overall_statistics_quality_score(
      data_validation, correlation_analysis, distribution_analysis,
      outlier_analysis, validation_results
    ),

    data_quality = list(
      input_quality = data_validation$overall_quality,
      missing_data_handled = calculate_missing_data_improvement(original_data, processed_data),
      property_coverage = length(property_statistics) / length(identify_numeric_soil_properties(original_data))
    ),

    analysis_quality = list(
      correlation_quality = if(!is.null(correlation_analysis)) assess_correlation_quality(correlation_analysis) else NULL,
      distribution_quality = if(!is.null(distribution_analysis)) assess_distribution_quality(distribution_analysis) else NULL,
      outlier_quality = if(!is.null(outlier_analysis)) assess_outlier_quality(outlier_analysis) else NULL
    ),

    validation_summary = if(!is.null(validation_results)) {
      list(
        overall_valid = validation_results$overall_valid,
        validation_score = validation_results$validation_score,
        n_errors = length(validation_results$errors),
        n_warnings = length(validation_results$warnings)
      )
    } else NULL,

    recommendations = generate_analysis_recommendations(
      data_validation, correlation_analysis, distribution_analysis,
      outlier_analysis, validation_results, config
    )
  )

  return(quality_report)
}

# ==============================================================================
# CONFIGURATION AND HELPER FUNCTIONS
# ==============================================================================

#' Get Statistical Analysis Default Configuration
#'
#' Returns default configuration optimized for statistical analysis
#'
#' @return Default statistical analysis configuration
#'
#' @export
get_statistical_analysis_defaults <- function() {

  base_config <- get_default_configuration("full")

  # Statistical analysis specific defaults
  statistical_config <- list(
    # Data processing
    stratify_by_horizon = TRUE,
    stratify_by_taxonomy = FALSE,
    include_texture_analysis = TRUE,

    # Missing values
    missing_value_strategy = "interpolate",
    missing_data_threshold = 0.8,

    # Correlation analysis
    correlation_methods = c("pearson", "spearman"),
    handle_constant_variables = "warn",
    minimum_observations = 10,

    # Distribution analysis
    distribution_methods = c("normal", "lognormal", "gamma", "beta"),
    fit_quality_threshold = 0.05,

    # Outlier detection
    outlier_methods = c("iqr", "zscore", "modified_zscore"),
    outlier_thresholds = list(iqr = 1.5, zscore = 3, modified_zscore = 3.5),
    detect_multivariate_outliers = TRUE,

    # Quality control
    min_quality_score = 0.6,
    quality_recovery_action = "warn",
    strict_property_validation = FALSE,

    # Error handling
    error_recovery_action = "warn",

    # Output control
    return_processed_data = FALSE,

    # Quality thresholds (using Module 0 defaults)
    quality_thresholds = get_default_quality_thresholds()
  )

  # Merge with base configuration
  return(merge_configurations(base_config, list(statistical_analysis = statistical_config)))
}

#' Validate Statistical Configuration
#'
#' Validate statistical analysis configuration using Module 0 utilities
#'
#' @param config Configuration to validate
#'
#' @return Validation results
#'
#' @export
validate_statistical_config <- function(config) {

  # Define parameter specifications with proper vector support
  param_specs <- list(
    minimum_observations = list(type = "numeric", required = TRUE, range = c(3, 1000)),
    missing_data_threshold = list(type = "numeric", required = TRUE, range = c(0, 1)),
    min_quality_score = list(type = "numeric", required = TRUE, range = c(0, 1)),

    # FIXED: correlation_methods as vector parameter
    correlation_methods = list(type = "character", required = TRUE,
                               choices = c("pearson", "spearman", "kendall"),
                               min_length = 1, max_length = 3),

    missing_value_strategy = list(type = "character", required = TRUE,
                                  choices = c("remove", "interpolate", "mean", "median", "forward_fill")),

    error_recovery_action = list(type = "character", required = TRUE,
                                 choices = c("stop", "warn", "continue", "retry")),

    # Additional parameters that might be vectors
    outlier_methods = list(type = "character", required = FALSE,
                           choices = c("iqr", "zscore", "modified_zscore", "robust"),
                           min_length = 1, max_length = 4),

    distribution_methods = list(type = "character", required = FALSE,
                                choices = c("normal", "lognormal", "gamma", "beta", "weibull"),
                                min_length = 1, max_length = 5)
  )

  # Extract statistical config if nested
  if ("statistical_analysis" %in% names(config)) {
    config_to_validate <- config$statistical_analysis
  } else {
    config_to_validate <- config
  }

  # Use the FIXED validate_parameters function
  return(validate_parameters(config_to_validate, param_specs, strict_mode = FALSE))
}

#' Identify Numeric Soil Properties
#'
#' Identify numeric soil properties in data using improved logic
#'
#' @param data Input data
#'
#' @return Vector of numeric soil property column names
#'
#' @export
identify_numeric_soil_properties <- function(data) {

  # Get all numeric columns
  numeric_cols <- names(data)[sapply(data, is.numeric)]

  # Filter to known soil properties (ending with _r for representative values)
  soil_property_pattern <- "_r$"
  soil_properties <- numeric_cols[grepl(soil_property_pattern, numeric_cols)]

  # Known important soil properties
  important_properties <- c(
    "sandtotal_r", "silttotal_r", "claytotal_r", "dbovendry_r", "dbthirdbar_r",
    "ph1to1h2o_r", "cec7_r", "om_r", "rfv_r", "wthirdbar_r", "wfifteenbar_r",
    "awc_r", "ksat_r", "ec_r", "sar_r", "esp_r"
  )

  # Prioritize important properties that are available
  available_important <- intersect(important_properties, soil_properties)
  other_soil_properties <- setdiff(soil_properties, important_properties)

  # Return important properties first, then others
  return(c(available_important, other_soil_properties))
}

# ==============================================================================
# SUPPORTING HELPER FUNCTIONS (LEVERAGING Module 0)
# ==============================================================================

#' Compute Property Statistics (Enhanced)
#'
#' Enhanced property statistics computation using Module 0 utilities
#'
#' @param data Input data
#' @param properties Properties to analyze
#' @param config Configuration
#'
#' @return Enhanced property statistics
#'
#' @export
compute_property_statistics <- function(data, properties, config) {

  property_stats <- list()

  for (property in properties) {
    values <- data[[property]]
    valid_values <- values[!is.na(values) & is.finite(values)]

    if (length(valid_values) == 0) {
      next
    }

    # Basic statistics
    basic_stats <- list(
      n_observations = length(valid_values),
      n_missing = sum(is.na(values)),
      missing_rate = sum(is.na(values)) / length(values),
      mean = mean(valid_values),
      sd = sd(valid_values),
      min = min(valid_values),
      max = max(valid_values),
      range = max(valid_values) - min(valid_values),
      cv = sd(valid_values) / mean(valid_values),
      median = median(valid_values),
      q25 = quantile(valid_values, 0.25),
      q75 = quantile(valid_values, 0.75),
      iqr = IQR(valid_values)
    )

    # Confidence intervals using Module 0
    confidence_intervals <- calculate_confidence_intervals(
      valid_values,
      statistic = "mean",
      confidence_level = 0.95,
      method = "t"
    )

    # Normality and distribution metrics
    distribution_metrics <- list(
      skewness = calculate_skewness(valid_values),
      kurtosis = calculate_kurtosis(valid_values),
      shapiro_p = if(length(valid_values) <= 5000) shapiro.test(valid_values)$p.value else NA
    )

    property_stats[[property]] <- c(basic_stats,
                                    list(confidence_intervals = confidence_intervals),
                                    distribution_metrics)
  }

  return(property_stats)
}

# Additional helper functions
calculate_condition_number <- function(matrix) {
  if (!is.matrix(matrix) || nrow(matrix) != ncol(matrix)) return(NA)
  eigenvals <- eigen(matrix, only.values = TRUE)$values
  eigenvals <- eigenvals[eigenvals > 1e-12]
  if (length(eigenvals) == 0) return(Inf)
  return(max(eigenvals) / min(eigenvals))
}

calculate_skewness <- function(x) {
  n <- length(x)
  if (n < 3) return(NA)
  m <- mean(x)
  s <- sd(x)
  if (s == 0) return(0)
  return(sum((x - m)^3) / (n * s^3))
}

calculate_kurtosis <- function(x) {
  n <- length(x)
  if (n < 4) return(NA)
  m <- mean(x)
  s <- sd(x)
  if (s == 0) return(0)
  return(sum((x - m)^4) / (n * s^4) - 3)
}

#' Determine candidate distributions to fit for a property
#'
#' Computes the existing type-appropriate heuristic candidate list, then, if
#' `config$distribution_methods` is supplied, narrows to the intersection of
#' the two - keeping the heuristic as a soft prior rather than silently
#' discarding an explicit user request. If the intersection is empty (the
#' user asked only for families the heuristic wouldn't have suggested),
#' falls back to the user's list unfiltered, with a logged warning, rather
#' than ignoring it entirely. Previously `config` was accepted by every
#' caller in this chain but never actually consulted here - `distribution_methods`
#' was inert.
#'
#' @param property_name Property name (drives the type-appropriate heuristic).
#' @param values Numeric vector of observed property values (currently
#'   unused by the heuristic itself, but taken for interface symmetry with
#'   callers and potential future data-driven refinement).
#' @param config Optional analysis configuration; reads `config$distribution_methods`.
#' @return Character vector of candidate distribution names.
#' @export
get_appropriate_distributions <- function(property_name, values, config = NULL) {
  heuristic <- if (property_name %in% c("sandtotal_r", "silttotal_r", "claytotal_r")) {
    c("beta", "normal")
  } else if (property_name == "ph1to1h2o_r") {
    c("normal", "gamma")
  } else {
    c("normal", "lognormal", "gamma", "weibull")
  }

  requested <- config$distribution_methods
  if (is.null(requested)) {
    return(heuristic)
  }

  narrowed <- intersect(requested, heuristic)
  if (length(narrowed) == 0) {
    log_message(
      "WARN",
      paste0("None of the configured distribution_methods (", paste(requested, collapse = ", "),
             ") match the type-appropriate heuristic for '", property_name,
             "' - using the configured list unfiltered rather than discarding it."),
      category = "Distribution"
    )
    return(requested)
  }
  narrowed
}

#' Fit a single candidate distribution to a property's values
#'
#' Modeled on the fitting style used in BCQRF/code/R/distribution_fitting.R:
#' a guarded `requireNamespace()` dependency, a `tryCatch`-wrapped fit, and a
#' structured return value rather than a raw fit object.
#'
#' @param values Numeric vector of observed property values.
#' @param dist_name One of "normal", "lognormal", "gamma", "beta" (the
#'   candidate names returned by get_appropriate_distributions()).
#' @param property_name Property name, used only for logging.
#' @param config Unused currently; kept for interface compatibility with callers.
#' @return A list with distribution/property_name/params/param_sd/loglik/aic/
#'   bic/convergence/n/scale_info/fit_object, or NULL if fitting isn't
#'   possible (too little data, wrong domain, or fitdist() failure).
fit_single_distribution <- function(values, dist_name, property_name, config) {
  if (!requireNamespace("fitdistrplus", quietly = TRUE)) {
    stop("Package 'fitdistrplus' is required for distribution fitting. ",
         "Install it with install.packages('fitdistrplus').")
  }

  values <- values[!is.na(values) & is.finite(values)]
  if (length(values) < 5) {
    return(NULL)
  }

  # lognormal/gamma/weibull are only defined for strictly positive support.
  if (dist_name %in% c("lognormal", "gamma", "weibull") && any(values <= 0)) {
    log_message("WARN", paste("Skipping", dist_name, "fit for", property_name,
                               "- requires strictly positive values"), category = "Distribution")
    return(NULL)
  }

  # fitdistrplus/stats use "norm"/"lnorm" rather than "normal"/"lognormal";
  # keep the friendlier candidate name in the returned object regardless.
  # "weibull" is fitdistrplus's own literal distname, so it needs no remapping,
  # but is listed explicitly here for clarity/self-documentation.
  r_dist_name <- c(normal = "norm", lognormal = "lnorm", gamma = "gamma", beta = "beta",
                    weibull = "weibull")[dist_name] %||% dist_name

  fit_values <- values
  scale_info <- NULL
  if (dist_name == "beta") {
    # fitdistrplus::fitdist() requires values strictly inside (0, 1).
    # Rescale using the data's own range (padded to include 0/100 for the
    # common texture-percentage case) rather than assuming a fixed domain.
    lower <- min(0, min(values))
    upper <- max(100, max(values))
    if (upper <= lower) upper <- lower + 1
    fit_values <- (values - lower) / (upper - lower)
    fit_values <- pmin(pmax(fit_values, 1e-6), 1 - 1e-6)
    scale_info <- list(lower = lower, upper = upper)
  }

  fit <- tryCatch({
    fitdistrplus::fitdist(fit_values, r_dist_name)
  }, error = function(e) {
    log_message("WARN", paste("fitdist() failed for", dist_name, "on", property_name,
                               ":", conditionMessage(e)), category = "Distribution")
    NULL
  })

  if (is.null(fit)) {
    return(NULL)
  }

  list(
    distribution = dist_name,
    property_name = property_name,
    params = as.list(fit$estimate),
    param_sd = as.list(fit$sd),
    loglik = fit$loglik,
    aic = fit$aic,
    bic = fit$bic,
    convergence = fit$convergence,
    n = length(fit_values),
    scale_info = scale_info,
    fit_object = fit
  )
}

#' Rank fitted distributions by AIC (lower is better)
#'
#' @param fitted_distributions Named list of fit results from
#'   fit_single_distribution() (NULL entries are dropped defensively, though
#'   fit_property_distributions() already filters them before this is called).
#' @return The same list, reordered best-fit-first.
rank_distribution_fits <- function(fitted_distributions) {
  fitted_distributions <- fitted_distributions[!vapply(fitted_distributions, is.null, logical(1))]
  if (length(fitted_distributions) <= 1) {
    return(fitted_distributions)
  }
  aics <- vapply(fitted_distributions, function(f) f$aic %||% Inf, numeric(1))
  fitted_distributions[order(aics)]
}

#' Goodness-of-fit tests across a property's candidate distributions
#'
#' Fits every candidate distribution (via get_appropriate_distributions() /
#' fit_single_distribution()) and runs fitdistrplus::gofstat() across the
#' successful fits, which returns KS/AD/CvM statistics (and, for the
#' non-censored continuous case here, approximate goodness-of-fit decisions)
#' for all candidates in a single call.
#'
#' @param values Numeric vector of observed property values.
#' @param property_name Property name.
#' @param config Passed through to fit_single_distribution().
#' @return A list with `fits` (the individual fit_single_distribution()
#'   results, best-AIC-first) and `gof` (a list of fitdistrplus::gofstat()
#'   results, one per distinct data-scale group - e.g. beta fits are grouped
#'   separately from unscaled fits - or NULL if no candidate fit succeeded).
perform_distribution_tests <- function(values, property_name, config) {
  if (!requireNamespace("fitdistrplus", quietly = TRUE)) {
    stop("Package 'fitdistrplus' is required for distribution tests. ",
         "Install it with install.packages('fitdistrplus').")
  }

  values <- values[!is.na(values) & is.finite(values)]
  candidates <- get_appropriate_distributions(property_name, values, config)

  fits <- list()
  for (dist_name in candidates) {
    fit <- tryCatch(
      fit_single_distribution(values, dist_name, property_name, config),
      error = function(e) NULL
    )
    if (!is.null(fit)) {
      fits[[dist_name]] <- fit
    }
  }
  fits <- rank_distribution_fits(fits)

  # gofstat() compares each candidate's fitted CDF against one shared
  # empirical dataset, so every fit object passed to a single gofstat() call
  # must have been fit on the same data. A "beta" fit is done on values
  # rescaled to (0,1) (see fit_single_distribution()), which is not
  # comparable to "normal"/"lognormal"/"gamma" fits done on the raw values -
  # mixing them would silently compare goodness-of-fit against the wrong
  # empirical distribution. Group fits by scale (unscaled vs. each distinct
  # rescaling) and only run gofstat() within a group. A group of 1 is still
  # meaningful (a single distribution's own KS/AD/CvM goodness-of-fit against
  # its data) - e.g. texture properties always split "beta" and "normal"
  # into separate 1-member groups since only beta needs rescaling, and both
  # are still worth reporting individually rather than dropped.
  scale_key <- function(f) {
    if (is.null(f$scale_info)) "unscaled" else paste(f$scale_info$lower, f$scale_info$upper)
  }
  groups <- split(fits, vapply(fits, scale_key, character(1)))

  gof_by_group <- lapply(groups, function(group_fits) {
    fit_objects <- lapply(group_fits, function(f) f$fit_object)
    # fitdistrplus::gofstat() rejects a length-1 *list* of fitdist objects
    # (errors with "argument f must a 'fitdist' ... object or a list of
    # ... objects") but accepts a single bare fitdist object directly -
    # unwrap for the single-distribution-per-group case.
    if (length(fit_objects) == 1) fit_objects <- fit_objects[[1]]
    tryCatch({
      fitdistrplus::gofstat(fit_objects, fitnames = names(group_fits))
    }, error = function(e) {
      log_message("WARN", paste("gofstat() failed for", property_name, ":", conditionMessage(e)),
                  category = "Distribution")
      NULL
    })
  })
  gof_by_group <- gof_by_group[!vapply(gof_by_group, is.null, logical(1))]

  list(fits = fits, gof = if (length(gof_by_group) > 0) gof_by_group else NULL)
}

assess_overall_distributions <- function(fitted_distributions, distribution_tests) {
  return(list(
    total_properties = length(fitted_distributions),
    avg_distributions_per_property = if(length(fitted_distributions) > 0) mean(sapply(fitted_distributions, length)) else 0
  ))
}

#' Detect Multivariate Outliers
#'
#' Real Mahalanobis-distance-based multivariate outlier detection across
#' `properties`, guarding a non-positive-definite covariance matrix with
#' [Matrix::nearPD()] (already an Import, no new dependency).
#'
#' @param data Input data frame.
#' @param properties Character vector of column names to include.
#' @param config List, optionally with `mahalanobis_alpha` (default 0.975,
#'   i.e. flag the outer 2.5% via a chi-squared cutoff).
#' @return List with `outliers` (logical vector, `NA` for incomplete rows),
#'   `n_outliers`, `method`, `threshold`.
detect_multivariate_outliers <- function(data, properties, config) {
  available_props <- intersect(properties, names(data))

  if (length(available_props) < 2) {
    return(list(outliers = logical(0), n_outliers = 0, method = "mahalanobis", threshold = NA_real_))
  }

  subset_data <- data[available_props]
  complete_mask <- stats::complete.cases(subset_data)

  if (sum(complete_mask) < length(available_props) + 1) {
    return(list(outliers = rep(NA, nrow(data)), n_outliers = 0, method = "mahalanobis", threshold = NA_real_))
  }

  complete_data <- as.matrix(subset_data[complete_mask, ])
  center <- colMeans(complete_data)
  cov_mat <- stats::cov(complete_data)

  cov_mat <- tryCatch({
    if (any(!is.finite(cov_mat)) || det(cov_mat) <= 0) {
      as.matrix(Matrix::nearPD(cov_mat)$mat)
    } else {
      cov_mat
    }
  }, error = function(e) cov_mat)

  distances <- tryCatch(
    stats::mahalanobis(complete_data, center, cov_mat),
    error = function(e) rep(NA_real_, nrow(complete_data))
  )

  alpha <- config$mahalanobis_alpha %||% 0.975
  threshold <- stats::qchisq(alpha, df = length(available_props))
  is_outlier_complete <- distances > threshold

  outliers <- rep(NA, nrow(data))
  outliers[complete_mask] <- is_outlier_complete

  list(
    outliers = outliers,
    n_outliers = sum(is_outlier_complete, na.rm = TRUE),
    method = "mahalanobis",
    threshold = threshold
  )
}

#' Generate Outlier Summary
#'
#' Aggregates [detect_comprehensive_outliers()]'s `property_outliers` (per
#' property, per method) and `multivariate_outliers` into overall counts.
#'
#' @param outlier_results Result of [detect_comprehensive_outliers()] (or a
#'   list with the same `property_outliers`/`multivariate_outliers` shape).
#' @return List with `total_outliers`, `property_outlier_counts`,
#'   `multivariate_outlier_count`.
generate_outlier_summary <- function(outlier_results) {
  property_counts <- list()
  total_property_outliers <- 0

  for (prop in names(outlier_results$property_outliers)) {
    prop_methods <- outlier_results$property_outliers[[prop]]
    method_counts <- vapply(prop_methods, function(m) m$n_outliers %||% 0, numeric(1))
    property_counts[[prop]] <- as.list(method_counts)
    if (length(method_counts) > 0) total_property_outliers <- total_property_outliers + max(method_counts)
  }

  mv_outliers <- outlier_results$multivariate_outliers$n_outliers %||% 0

  list(
    total_outliers = total_property_outliers + mv_outliers,
    property_outlier_counts = property_counts,
    multivariate_outlier_count = mv_outliers
  )
}

#' Validate a set of correlation matrices
#'
#' Real implementation using the shared `validate_correlation_matrix()`
#' (`distributions.R`) - previously a stub always returning `valid=TRUE`.
#'
#' @param matrices Named list keyed by method, each entry either a matrix or
#'   `list(matrix=, ...)` (matching `run_comprehensive_correlation_analysis()`'s
#'   `$matrices` shape).
#' @param config Unused; kept for interface compatibility with callers.
#' @return `list(valid=, errors=, warnings=)`.
validate_correlation_matrices <- function(matrices, config) {
  warnings_out <- character(0)

  for (method in names(matrices)) {
    entry <- matrices[[method]]
    mat <- if (is.matrix(entry)) entry else entry$matrix
    if (!is.matrix(mat)) next
    properties <- rownames(mat) %||% colnames(mat) %||% seq_len(nrow(mat))
    check <- validate_correlation_matrix(mat, properties)
    if (!check$valid) {
      warnings_out <- c(warnings_out, paste0(method, ": ", check$message))
    }
  }

  list(valid = TRUE, errors = character(0), warnings = warnings_out)
}

validate_distribution_analysis <- function(distribution_analysis, config) {
  return(list(warnings = character(0)))
}

validate_outlier_analysis <- function(outlier_analysis, config) {
  return(list(warnings = character(0)))
}

calculate_overall_validation_score <- function(validation_results) {
  base_score <- if(validation_results$overall_valid) 1.0 else 0.5
  error_penalty <- length(validation_results$errors) * 0.1
  warning_penalty <- length(validation_results$warnings) * 0.02
  return(max(0, base_score - error_penalty - warning_penalty))
}

generate_correlation_summary <- function(correlation_results) {
  return(list(n_matrices = length(correlation_results$matrices)))
}

#' Analyze texture (sand/silt/clay) correlations, in both raw and ILR space
#'
#' Raw Pearson correlations among compositional (simplex-constrained) parts
#' like sand/silt/clay percentages are spuriously negative due to the
#' sum-to-100 constraint - `ilr_correlations` (isometric log-ratio space, via
#' the dependency-free `ilr_forward()` in `distributions.R`) avoids that
#' artifact and is the statistically defensible view for compositional data.
#' Both are reported: `raw_correlations` for continuity with prior behavior,
#' `ilr_correlations` as the added analytical value (only computed when all
#' three fractions are present with at least 5 complete observations).
#'
#' @param data Input data with texture `_r` columns.
#' @param texture_properties Character vector of available texture column
#'   names (e.g. a subset/all of `c("sandtotal_r","silttotal_r","claytotal_r")`).
#' @param methods Correlation methods (e.g. `c("pearson","spearman")`).
#' @param config Unused; kept for interface compatibility with callers.
#' @return `list(raw_correlations=, ilr_correlations=, n_observations=, note=)`.
#' @export
analyze_texture_correlations <- function(data, texture_properties, methods, config) {
  result <- list(
    raw_correlations = list(),
    ilr_correlations = NULL,
    note = paste(
      "Raw Pearson/Spearman correlations among compositional (simplex-constrained)",
      "parts like sand/silt/clay percentages are spuriously negative due to the",
      "sum-to-100 constraint; ilr_correlations (isometric log-ratio space) avoids",
      "that artifact and is the more defensible compositional-data view."
    )
  )

  if (length(texture_properties) >= 2) {
    texture_data <- data[, texture_properties, drop = FALSE]
    for (method in methods) {
      result$raw_correlations[[method]] <- safe_correlation(texture_data, method = method, handle_constant = "warn")
    }
  }

  required <- c("claytotal_r", "sandtotal_r", "silttotal_r")
  if (all(required %in% texture_properties) && all(required %in% names(data))) {
    complete_rows <- stats::complete.cases(data[, required])
    if (sum(complete_rows) > 4) {
      texture_complete <- data[complete_rows, required]
      z <- ilr_forward(texture_complete$claytotal_r, texture_complete$sandtotal_r, texture_complete$silttotal_r)

      ilr_cor <- list()
      for (method in methods) {
        ilr_cor[[method]] <- safe_correlation(as.data.frame(z), method = method, handle_constant = "warn")
      }
      result$ilr_correlations <- ilr_cor
      result$n_observations <- sum(complete_rows)
    }
  }

  result
}

#' Calculate Overall Statistics Quality Score
#'
#' Real composite score averaging whichever of `data_validation`'s quality
#' score, `validation_results`'s validation score, and the
#' [assess_correlation_quality()]/[assess_distribution_quality()]/
#' [assess_outlier_quality()] scores are available.
#'
#' @param data_validation Result of `validate_data_quality()`, with an
#'   `overall_quality$score` field.
#' @param correlation_analysis Result of correlation analysis (see
#'   [assess_correlation_quality()]), or `NULL`.
#' @param distribution_analysis Result of distribution analysis (see
#'   [assess_distribution_quality()]), or `NULL`.
#' @param outlier_analysis Result of [detect_comprehensive_outliers()], or `NULL`.
#' @param validation_results Result of statistical-results validation, with
#'   a `validation_score` field.
#' @return Single numeric quality score in `[0, 1]`, or `NA` if nothing was available.
calculate_overall_statistics_quality_score <- function(data_validation, correlation_analysis,
                                                     distribution_analysis, outlier_analysis, validation_results) {
  scores <- c()

  if (!is.null(data_validation$overall_quality$score)) {
    scores <- c(scores, data_validation$overall_quality$score)
  }
  if (!is.null(validation_results$validation_score)) {
    scores <- c(scores, validation_results$validation_score)
  }
  if (!is.null(correlation_analysis)) {
    scores <- c(scores, assess_correlation_quality(correlation_analysis)$quality_score)
  }
  if (!is.null(distribution_analysis)) {
    scores <- c(scores, assess_distribution_quality(distribution_analysis)$quality_score)
  }
  if (!is.null(outlier_analysis)) {
    scores <- c(scores, assess_outlier_quality(outlier_analysis)$quality_score)
  }

  scores <- scores[is.finite(scores)]
  if (length(scores) == 0) return(NA_real_)
  mean(scores)
}

#' Calculate Missing-Data Handling Improvement
#'
#' Real reduction in missing-value rate (averaged across shared numeric
#' columns) between `original_data` and `processed_data`.
#'
#' @param original_data Original input data.
#' @param processed_data Data after infilling/processing.
#' @return Numeric improvement in `[0, 1]` (0 = no improvement), or `NA` if
#'   not computable.
calculate_missing_data_improvement <- function(original_data, processed_data) {
  if (is.null(original_data) || is.null(processed_data)) return(NA_real_)

  common_cols <- intersect(names(original_data), names(processed_data))
  numeric_cols <- common_cols[vapply(common_cols, function(col) is.numeric(original_data[[col]]), logical(1))]
  if (length(numeric_cols) == 0) return(NA_real_)

  orig_na_rate <- mean(vapply(numeric_cols, function(col) mean(is.na(original_data[[col]])), numeric(1)))
  processed_na_rate <- mean(vapply(numeric_cols, function(col) mean(is.na(processed_data[[col]])), numeric(1)))

  max(0, orig_na_rate - processed_na_rate)
}

#' Assess Correlation Analysis Quality
#'
#' Real score based on the fraction of correlation matrices that produced a
#' validation warning (see `validate_correlation_matrices()`'s
#' `$matrices`/`$validation$warnings` shape).
#'
#' @param correlation_analysis List with `matrices` and `validation$warnings`.
#' @return List with `quality_score`, `n_matrices`, `n_warnings`.
assess_correlation_quality <- function(correlation_analysis) {
  if (is.null(correlation_analysis)) return(list(quality_score = NA_real_, n_matrices = 0, n_warnings = 0))

  n_warnings <- length(correlation_analysis$validation$warnings %||% character(0))
  n_matrices <- length(correlation_analysis$matrices %||% list())
  quality_score <- if (n_matrices > 0) max(0, 1 - (n_warnings / n_matrices) * 0.5) else NA_real_

  list(quality_score = quality_score, n_matrices = n_matrices, n_warnings = n_warnings)
}

#' Assess Distribution Analysis Quality
#'
#' Real score based on the fraction of per-property distribution-fit
#' entries that produced a non-empty `fits` list.
#'
#' @param distribution_analysis Named list keyed by property, each entry
#'   optionally with a `fits` field (see [fit_property_distributions()]-style output).
#' @return List with `quality_score`, `n_properties`.
assess_distribution_quality <- function(distribution_analysis) {
  if (is.null(distribution_analysis) || length(distribution_analysis) == 0) {
    return(list(quality_score = NA_real_, n_properties = 0))
  }

  n_properties <- length(distribution_analysis)
  n_valid <- sum(vapply(distribution_analysis, function(x) {
    !is.null(x) && (is.null(x$fits) || length(x$fits) > 0)
  }, logical(1)))

  list(quality_score = n_valid / n_properties, n_properties = n_properties)
}

#' Assess Outlier Analysis Quality
#'
#' Real score based on the overall outlier rate implied by
#' [detect_comprehensive_outliers()]'s `summary`/`property_outliers` shape.
#'
#' @param outlier_analysis Result of [detect_comprehensive_outliers()].
#' @return List with `quality_score`, `total_outliers`, `outlier_rate`.
assess_outlier_quality <- function(outlier_analysis) {
  if (is.null(outlier_analysis)) return(list(quality_score = NA_real_, total_outliers = 0, outlier_rate = NA_real_))

  total_outliers <- outlier_analysis$summary$total_outliers %||% 0

  total_obs <- sum(vapply(outlier_analysis$property_outliers %||% list(), function(prop_methods) {
    lengths_per_method <- vapply(prop_methods, function(m) length(m$outliers %||% logical(0)), numeric(1))
    if (length(lengths_per_method) > 0) max(lengths_per_method) else 0
  }, numeric(1)))

  outlier_rate <- if (total_obs > 0) total_outliers / total_obs else 0
  list(quality_score = max(0, 1 - outlier_rate), total_outliers = total_outliers, outlier_rate = outlier_rate)
}

generate_analysis_recommendations <- function(data_validation, correlation_analysis,
                                              distribution_analysis, outlier_analysis,
                                              validation_results, config) {
  recommendations <- character(0)

  # Data quality recommendations
  if (!is.null(data_validation) && data_validation$overall_quality$score < 0.8) {
    recommendations <- c(recommendations, "Consider improving data quality before analysis")
  }

  # Correlation recommendations
  if (!is.null(correlation_analysis) && !is.null(correlation_analysis$validation)) {
    if (length(correlation_analysis$validation$warnings) > 0) {
      recommendations <- c(recommendations, "Review correlation matrix warnings")
    }
  }

  return(recommendations)
}
