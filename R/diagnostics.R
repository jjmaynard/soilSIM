#' @title Workflow Validation & Diagnostics
#' @description Functions for validating and diagnosing the full soil
#'   simulation workflow - Monte Carlo output quality, correlation structure
#'   preservation, GP model performance, and soil-science plausibility -
#'   plus diagnostic plot generation and reporting.
#' @name validation_diagnostics
NULL

# ============================================================================
# 1. MASTER VALIDATION FUNCTIONS
# ============================================================================

#' Validate Complete Workflow
#'
#' Assesses the entire soil simulation workflow end to end.
#'
#' @param workflow_results Complete workflow results including all intermediate steps
#' @param original_data Original SSURGO/NRCS data for comparison
#' @param validation_config Configuration for validation parameters (uses package defaults)
#' @param generate_plots Whether to generate diagnostic plots (default = TRUE)
#' @param output_dir Directory for saving validation outputs
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's log level so \code{INFO}-level progress messages print for the duration of this call (default \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @return Comprehensive validation results, as a `soilSIM_diagnostics` object
#' @export
diagnose_workflow <- function(workflow_results,
                                       original_data = NULL,
                                       validation_config = NULL,
                                       generate_plots = TRUE,
                                       output_dir = NULL,
                                       verbose = getOption("ssurgo.verbose", FALSE)) {

  .old_log_cfg <- set_verbose_logging(verbose)
  on.exit(options(soil_workflow_log_config = .old_log_cfg), add = TRUE)

  log_message("INFO", "=== COMPREHENSIVE WORKFLOW VALIDATION ===", category = "Validation")
  start_time <- Sys.time()

  # Use the shared configuration helpers
  if (is.null(validation_config)) {
    validation_config <- default_config("validation")
  }

  # Initialize validation results structure
  validation_results <- initialize_validation_structure(start_time, validation_config, output_dir)

  # Extract and validate workflow components utilities
  components <- extract_and_validate_components(workflow_results)

  log_message("INFO", paste("Validating workflow with", length(components), "components"), category = "Validation")

  # validation pipeline
  validation_results <- execute_validation_pipeline(
    validation_results, components, original_data, validation_config, generate_plots, output_dir
  )

  # Finalize validation
  validation_results <- finalize_validation_results(validation_results, start_time)

  log_message("INFO", "=== VALIDATION COMPLETE ===", category = "Validation")
  log_message("INFO", paste("Overall quality score:", round(validation_results$overall_assessment$quality_score, 3)), category = "Validation")

  return(new_soilSIM_diagnostics(validation_results))
}

#' Generate Validation Report
#'
#'
#' @param validation_results Results from diagnose_workflow()
#' @param output_format Format for report: "html", "pdf", or "markdown"
#' @param output_file Output file path
#' @param include_plots Whether to include diagnostic plots
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's log level so \code{INFO}-level progress messages print for the duration of this call (default \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @return Success status
#' @keywords internal
generate_validation_report <- function(validation_results,
                                       output_format = "html",
                                       output_file = NULL,
                                       include_plots = TRUE,
                                       verbose = getOption("ssurgo.verbose", FALSE)) {

  .old_log_cfg <- set_verbose_logging(verbose)
  on.exit(options(soil_workflow_log_config = .old_log_cfg), add = TRUE)

  log_message("INFO", "=== GENERATING VALIDATION REPORT ===", category = "Validation")
  log_message("INFO", paste("Output format:", output_format), category = "Validation")

  # file naming utilities
  if (is.null(output_file)) {
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    output_file <- paste0("soil_simulation_validation_", timestamp, ".", output_format)
  }

  log_message("INFO", paste("Output file:", output_file), category = "Validation")

  # Validate parameters
  param_specs <- list(
    output_format = list(required = TRUE, type = "character",
                         choices = c("html", "pdf", "markdown")),
    include_plots = list(required = FALSE, type = "logical")
  )

  param_validation <- validate_parameters(
    list(output_format = output_format, include_plots = include_plots),
    param_specs,
    strict_mode = TRUE
  )

  if (!param_validation$valid) {
    log_message("ERROR", paste("Parameter validation failed:", paste(param_validation$errors, collapse = ", ")), category = "Validation")
    return(FALSE)
  }

  # Create report content with error handling
  report_content <- tryCatch({
    create_report_content(validation_results, include_plots)
  }, error = function(e) {
    handle_workflow_error(e, "Report content creation", "warn")
    return(NULL)
  })

  if (is.null(report_content)) {
    log_message("ERROR", "Failed to create report content", category = "Validation")
    return(FALSE)
  }

  # Generate report
  success <- tryCatch({
    switch(output_format,
           "html" = generate_html_report(report_content, output_file),
           "pdf" = generate_pdf_report(report_content, output_file),
           "markdown" = generate_markdown_report(report_content, output_file),
           {
             log_message("ERROR", paste("Unsupported output format:", output_format), category = "Validation")
             FALSE
           }
    )
  }, error = function(e) {
    handle_workflow_error(e, paste("Report generation for", output_format), "warn")
    return(FALSE)
  })

  if (success) {
    log_message("INFO", "Validation report generated successfully", category = "Validation")
  } else {
    log_message("ERROR", "Report generation failed", category = "Validation")
  }

  return(success)
}

#' Assess Workflow Quality
#'
#'
#' @param validation_results Complete validation results
#' @param validation_config Validation configuration
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's log level so \code{INFO}-level progress messages print for the duration of this call (default \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @return Overall quality assessment
#' @keywords internal
assess_workflow_quality <- function(validation_results, validation_config, verbose = getOption("ssurgo.verbose", FALSE)) {

  .old_log_cfg <- set_verbose_logging(verbose)
  on.exit(options(soil_workflow_log_config = .old_log_cfg), add = TRUE)

  log_message("INFO", "Performing overall workflow quality assessment", category = "Validation")

  # Component quality scores safe operations
  quality_scores <- calculate_component_quality_scores(validation_results)

  # Calculate weighted overall score utilities
  weights <- validation_config$quality_weights %||% get_default_quality_weights()
  overall_score <- calculate_weighted_quality_score(quality_scores, weights)

  # Determine quality grade patterns
  quality_grade <- determine_quality_grade(overall_score)

  # Identify critical issues validation framework
  critical_issues <- identify_critical_issues(validation_results, validation_config)

  # Create assessment
  assessment <- create_quality_assessment(
    overall_score, quality_grade, quality_scores, critical_issues, validation_config
  )

  log_message("INFO", paste("Quality assessment - Score:", round(overall_score, 3),
                            "Grade:", quality_grade, "Status:", assessment$workflow_status), category = "Validation")

  return(assessment)
}

# ============================================================================
# 2. MONTE CARLO VALIDATION FUNCTIONS# ============================================================================

#' Validate Monte Carlo Quality
#'
#'
#' @param monte_carlo_results Results from monte_carlo module
#' @param original_data Original SSURGO data for comparison
#' @param config Monte Carlo validation configuration
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's log level so \code{INFO}-level progress messages print for the duration of this call (default \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @return Monte Carlo validation results
#' @keywords internal
diagnose_simulation <- function(monte_carlo_results,
                                         original_data = NULL,
                                         config = NULL,
                                         verbose = getOption("ssurgo.verbose", FALSE)) {

  .old_log_cfg <- set_verbose_logging(verbose)
  on.exit(options(soil_workflow_log_config = .old_log_cfg), add = TRUE)

  log_message("INFO", "Validating Monte Carlo simulation quality", category = "Validation")

  # Use the shared configuration helpers
  if (is.null(config)) {
    config <- default_config("monte_carlo")
  }

  # data validation
  simulation_data <- monte_carlo_results$simulation_data

  data_validation <- validate_data_quality(
    simulation_data,
    required_columns = c("cokey", "simulation_number"),
    quality_thresholds = config$data_quality_thresholds
  )

  if (!data_validation$overall_quality$validation_passed) {
    log_message("WARN", "Monte Carlo simulation data has quality issues", category = "Validation")
  }

  # Initialize validation structure
  validation <- list(
    convergence_assessment = list(),
    coverage_assessment = list(),
    distribution_fidelity = list(),
    simulation_diagnostics = list(),
    data_quality = data_validation
  )

  # validation pipeline with progress tracking
  validation_steps <- list(
    "convergence" = function() assess_simulation_convergence(simulation_data, config$convergence_criteria),
    "coverage" = function() diagnose_simulation_coverage(simulation_data, original_data, config$coverage_criteria),
    "distribution" = function() diagnose_distribution_fidelity(simulation_data, monte_carlo_results$metadata, config$distribution_criteria),
    # generate_simulation_diagnostics() takes 5 args (original_data, simulation_results,
    # properties, correlation_config, config) - a prior 2-arg call here always errored and was
    # silently swallowed by the tryCatch below into a validation_failed placeholder. Sourced from
    # monte_carlo_results, mirroring simulate_monte_carlo()'s own call convention.
    "diagnostics" = function() generate_simulation_diagnostics(
      original_data,
      simulation_data,
      monte_carlo_results$metadata$properties,
      monte_carlo_results$correlation_structure,
      config
    )
  )

  for (i in seq_along(validation_steps)) {
    step_name <- names(validation_steps)[i]
    log_message("DEBUG", paste("Assessing", step_name), category = "Validation")

    track_progress(i, length(validation_steps), "Monte Carlo validation", update_frequency = 1)

    validation[[paste0(step_name, "_assessment")]] <- tryCatch({
      validation_steps[[i]]()
    }, error = function(e) {
      handle_workflow_error(e, paste("Monte Carlo", step_name, "validation"), "warn")
      return(list(validation_failed = TRUE, error = e$message))
    })
  }

  log_message("INFO", "Monte Carlo validation completed", category = "Validation")
  return(validation)
}

#' Assess Simulation Coverage
#'
#'
#' @param simulation_data Simulation results
#' @param original_data Original data for comparison
#' @param criteria Coverage assessment criteria
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's log level so \code{INFO}-level progress messages print for the duration of this call (default \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @return Coverage assessment results
#' @keywords internal
diagnose_simulation_coverage <- function(simulation_data,
                                                original_data,
                                                criteria = NULL,
                                                verbose = getOption("ssurgo.verbose", FALSE)) {

  .old_log_cfg <- set_verbose_logging(verbose)
  on.exit(options(soil_workflow_log_config = .old_log_cfg), add = TRUE)

  if (is.null(criteria)) {
    criteria <- list(
      min_coverage_percentile = 0.95,
      max_extrapolation_factor = 1.2,
      min_samples_per_group = 10
    )
  }

  log_message("DEBUG", "Assessing simulation coverage", category = "Validation")

  coverage_results <- list(
    parameter_space_coverage = list(),
    distributional_coverage = list(),
    group_representation = list()
  )

  # property detection
  properties_validation <- validate_properties(names(simulation_data), "laboratory", strict_mode = FALSE)
  numeric_properties <- intersect(
    properties_validation$property_stats$valid_property_names %||% names(simulation_data),
    names(simulation_data)[sapply(simulation_data, is.numeric)]
  )

  # Remove metadata columns
  numeric_properties <- setdiff(numeric_properties, c("simulation_number", "hzdept_r", "hzdepb_r"))

  log_message("DEBUG", paste("Assessing coverage for", length(numeric_properties), "properties"), category = "Validation")

  # Parameter space coverage with error handling
  if (!is.null(original_data)) {
    for (prop in numeric_properties) {
      if (prop %in% names(original_data)) {
        coverage_results$parameter_space_coverage[[prop]] <- tryCatch({
          assess_property_coverage(
            simulation_data[[prop]],
            original_data[[prop]],
            criteria
          )
        }, error = function(e) {
          handle_workflow_error(e, paste("Coverage assessment for", prop), "warn")
          return(list(coverage_failed = TRUE))
        })
      }
    }
  }

  # Distributional coverage safe operations
  for (prop in numeric_properties) {
    coverage_results$distributional_coverage[[prop]] <- tryCatch({
      assess_distributional_coverage(simulation_data[[prop]], criteria)
    }, error = function(e) {
      handle_workflow_error(e, paste("Distributional coverage for", prop), "warn")
      return(list(coverage_failed = TRUE))
    })
  }

  # Group representation
  if ("cokey" %in% names(simulation_data)) {
    coverage_results$group_representation <- tryCatch({
      assess_group_representation(simulation_data, criteria)
    }, error = function(e) {
      handle_workflow_error(e, "Group representation assessment", "warn")
      return(list(assessment_failed = TRUE))
    })
  }

  return(coverage_results)
}

#' Validate Distribution Fidelity
#'
#'
#' @param simulation_data Simulation results
#' @param simulation_metadata Metadata from Monte Carlo generation
#' @param criteria Distribution validation criteria
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's log level so \code{INFO}-level progress messages print for the duration of this call (default \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @return Distribution fidelity assessment
#' @keywords internal
diagnose_distribution_fidelity <- function(simulation_data,
                                                    simulation_metadata,
                                                    criteria = NULL,
                                                    verbose = getOption("ssurgo.verbose", FALSE)) {

  .old_log_cfg <- set_verbose_logging(verbose)
  on.exit(options(soil_workflow_log_config = .old_log_cfg), add = TRUE)

  if (is.null(criteria)) {
    criteria <- list(
      ks_test_alpha = 0.05,
      moment_tolerance = 0.1,
      quantile_tolerance = 0.05
    )
  }

  log_message("DEBUG", "Validating distribution fidelity", category = "Validation")

  fidelity_results <- list(
    distribution_tests = list(),
    moment_comparisons = list(),
    quantile_comparisons = list(),
    outlier_assessment = list()
  )

  # distribution parameter extraction
  if (!is.null(simulation_metadata$distribution_parameters)) {

    for (prop in names(simulation_metadata$distribution_parameters)) {
      if (prop %in% names(simulation_data)) {

        param_info <- simulation_metadata$distribution_parameters[[prop]]
        simulated_values <- simulation_data[[prop]][!is.na(simulation_data[[prop]])]

        if (length(simulated_values) > 0) {

          # Distribution tests
          fidelity_results$distribution_tests[[prop]] <- tryCatch({
            test_distribution_fidelity(simulated_values, param_info, criteria)
          }, error = function(e) {
            handle_workflow_error(e, paste("Distribution test for", prop), "warn")
            return(list(test_failed = TRUE))
          })

          # Moment comparisons safe operations
          fidelity_results$moment_comparisons[[prop]] <- tryCatch({
            compare_distribution_moments(simulated_values, param_info, criteria)
          }, error = function(e) {
            handle_workflow_error(e, paste("Moment comparison for", prop), "warn")
            return(list(comparison_failed = TRUE))
          })

          # Quantile comparisons
          fidelity_results$quantile_comparisons[[prop]] <- tryCatch({
            compare_distribution_quantiles(simulated_values, param_info, criteria)
          }, error = function(e) {
            handle_workflow_error(e, paste("Quantile comparison for", prop), "warn")
            return(list(comparison_failed = TRUE))
          })

          # Outlier assessment utilities
          fidelity_results$outlier_assessment[[prop]] <- tryCatch({
            outliers <- detect_outliers(simulated_values, method = "iqr", return_indices = FALSE)
            list(
              n_outliers = sum(outliers, na.rm = TRUE),
              outlier_percentage = mean(outliers, na.rm = TRUE) * 100,
              outlier_threshold_exceeded = mean(outliers, na.rm = TRUE) > criteria$outlier_threshold %||% 0.05
            )
          }, error = function(e) {
            handle_workflow_error(e, paste("Outlier assessment for", prop), "warn")
            return(list(assessment_failed = TRUE))
          })
        }
      }
    }
  }

  return(fidelity_results)
}

# ============================================================================
# 3. CORRELATION STRUCTURE VALIDATION
# ============================================================================

#' Validate Correlation Structures
#'
#'
#' @param correlation_matrices Correlation matrices from correlation_structure module
#' @param simulation_data Final simulation data
#' @param config Correlation validation configuration
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's log level so \code{INFO}-level progress messages print for the duration of this call (default \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @return Correlation validation results
#' @keywords internal
diagnose_correlation_structures <- function(correlation_matrices,
                                            simulation_data,
                                            config = NULL,
                                            verbose = getOption("ssurgo.verbose", FALSE)) {

  .old_log_cfg <- set_verbose_logging(verbose)
  on.exit(options(soil_workflow_log_config = .old_log_cfg), add = TRUE)

  log_message("INFO", "Validating correlation structures", category = "Validation")

  # Use the shared configuration helpers
  if (is.null(config)) {
    config <- default_config("correlation")
  }

  # data validation
  correlation_data_validation <- validate_data_quality(
    simulation_data,
    required_columns = c("cokey", "simulation_number"),
    numeric_columns = names(simulation_data)[sapply(simulation_data, is.numeric)]
  )

  validation <- list(
    matrix_quality = list(),
    preservation_assessment = list(),
    depth_specific_validation = list(),
    cholesky_validation = list(),
    data_quality = correlation_data_validation
  )

  # validation pipeline
  validation_steps <- list(
    "matrix_quality" = function() validate_correlation_matrix_quality(correlation_matrices, config$matrix_criteria),
    "preservation" = function() diagnose_correlation_preservation(correlation_matrices, simulation_data, config$preservation_criteria),
    "depth_specific" = function() diagnose_within_depth_correlations(simulation_data, config$depth_criteria),
    "cholesky" = function() diagnose_cholesky(correlation_matrices, config$cholesky_criteria)
  )

  for (i in seq_along(validation_steps)) {
    step_name <- names(validation_steps)[i]
    log_message("DEBUG", paste("Validating", step_name), category = "Validation")

    track_progress(i, length(validation_steps), "Correlation validation", update_frequency = 1)

    validation[[step_name]] <- tryCatch({
      validation_steps[[i]]()
    }, error = function(e) {
      handle_workflow_error(e, paste("Correlation", step_name, "validation"), "warn")
      return(list(validation_failed = TRUE, error = e$message))
    })
  }

  log_message("INFO", "Correlation structure validation completed", category = "Validation")
  return(validation)
}

#' Validate Correlation Preservation
#'
#'
#' @param original_correlations Original correlation matrices
#' @param simulation_data Final simulation data
#' @param criteria Preservation criteria
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's log level so \code{INFO}-level progress messages print for the duration of this call (default \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @return Correlation preservation assessment
#' @keywords internal
diagnose_correlation_preservation <- function(original_correlations,
                                               simulation_data,
                                               criteria = NULL,
                                               verbose = getOption("ssurgo.verbose", FALSE)) {

  .old_log_cfg <- set_verbose_logging(verbose)
  on.exit(options(soil_workflow_log_config = .old_log_cfg), add = TRUE)

  if (is.null(criteria)) {
    criteria <- list(
      max_correlation_difference = 0.1,
      correlation_rmse_threshold = 0.05,
      min_samples_for_validation = 30
    )
  }

  log_message("DEBUG", "Assessing correlation preservation", category = "Validation")

  preservation_results <- list(
    overall_preservation = list(),
    depth_specific_preservation = list(),
    property_specific_preservation = list()
  )

  # property detection
  properties_validation <- validate_properties(names(simulation_data), "laboratory", strict_mode = FALSE)
  numeric_properties <- intersect(
    properties_validation$property_stats$valid_property_names %||% names(simulation_data),
    names(simulation_data)[sapply(simulation_data, is.numeric)]
  )

  numeric_properties <- setdiff(numeric_properties, c("simulation_number", "hzdept_r", "hzdepb_r"))

  if (length(numeric_properties) < 2) {
    log_message("WARN", "Insufficient numeric properties for correlation assessment", category = "Validation")
    return(preservation_results)
  }

  # Overall preservation assessment safe correlation
  overall_sim_cor <- tryCatch({
    safe_correlation(simulation_data[numeric_properties], method = "pearson", handle_constant = "warn")
  }, error = function(e) {
    handle_workflow_error(e, "Overall correlation calculation", "warn")
    return(NULL)
  })

  if (!is.null(overall_sim_cor) && !is.null(original_correlations$global_correlation_matrix)) {
    original_cor <- original_correlations$global_correlation_matrix

    # Match dimensions safely
    common_props <- intersect(rownames(original_cor), rownames(overall_sim_cor))

    if (length(common_props) >= 2) {
      orig_subset <- original_cor[common_props, common_props]
      sim_subset <- overall_sim_cor[common_props, common_props]

      preservation_results$overall_preservation <- assess_correlation_differences(
        orig_subset, sim_subset, criteria
      )
    }
  }

  # Depth-specific preservation with error handling
  unique_depths <- unique(simulation_data$hzdept_r)
  depths_to_check <- unique_depths[1:min(5, length(unique_depths))]

  for (depth in depths_to_check) {
    depth_data <- simulation_data |> dplyr::filter(hzdept_r == depth)

    if (nrow(depth_data) >= criteria$min_samples_for_validation) {
      depth_cor <- tryCatch({
        safe_correlation(depth_data[numeric_properties], method = "pearson", handle_constant = "warn")
      }, error = function(e) {
        handle_workflow_error(e, paste("Depth correlation calculation for depth", depth), "warn")
        return(NULL)
      })

      if (!is.null(depth_cor)) {
        preservation_results$depth_specific_preservation[[as.character(depth)]] <-
          assess_depth_correlation_preservation(depth_cor, original_correlations, depth, criteria)
      }
    }
  }

  return(preservation_results)
}

#' Validate Within Depth Correlations
#'
#'
#' @param simulation_data Simulation data
#' @param criteria Depth validation criteria
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's log level so \code{INFO}-level progress messages print for the duration of this call (default \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @return Within-depth correlation validation
#' @keywords internal
diagnose_within_depth_correlations <- function(simulation_data, criteria = NULL, verbose = getOption("ssurgo.verbose", FALSE)) {

  .old_log_cfg <- set_verbose_logging(verbose)
  on.exit(options(soil_workflow_log_config = .old_log_cfg), add = TRUE)

  if (is.null(criteria)) {
    criteria <- list(
      depth_bins = c(0, 15, 30, 60, 100, 200),
      min_samples_per_bin = 20,
      expected_correlation_ranges = list()
    )
  }

  log_message("DEBUG", "Validating within-depth correlations", category = "Validation")

  depth_validation <- list(
    depth_bin_correlations = list(),
    correlation_stability = list(),
    depth_trend_correlations = list()
  )

  # property detection
  properties_validation <- validate_properties(names(simulation_data), "laboratory", strict_mode = FALSE)
  numeric_properties <- intersect(
    properties_validation$property_stats$valid_property_names %||% names(simulation_data),
    names(simulation_data)[sapply(simulation_data, is.numeric)]
  )

  numeric_properties <- setdiff(numeric_properties, c("simulation_number", "hzdept_r", "hzdepb_r"))

  if (length(numeric_properties) < 2) {
    log_message("WARN", "Insufficient numeric properties for within-depth correlation validation", category = "Validation")
    return(depth_validation)
  }

  # Assess correlations within depth bins safe operations
  simulation_data$depth_bin <- cut(simulation_data$hzdept_r,
                                   breaks = criteria$depth_bins,
                                   include.lowest = TRUE)

  for (bin in levels(simulation_data$depth_bin)) {
    bin_data <- simulation_data |> dplyr::filter(depth_bin == bin, !is.na(depth_bin))

    if (nrow(bin_data) >= criteria$min_samples_per_bin) {
      bin_cor <- tryCatch({
        safe_correlation(bin_data[numeric_properties], method = "pearson", handle_constant = "warn")
      }, error = function(e) {
        handle_workflow_error(e, paste("Bin correlation calculation for", bin), "warn")
        return(NULL)
      })

      if (!is.null(bin_cor)) {
        depth_validation$depth_bin_correlations[[bin]] <- list(
          correlation_matrix = bin_cor,
          n_observations = nrow(bin_data),
          mean_absolute_correlation = mean(abs(bin_cor[upper.tri(bin_cor)]), na.rm = TRUE),
          correlation_quality = assess_correlation_matrix_properties(bin_cor)
        )
      }
    }
  }

  # Assess correlation stability across depths
  depth_validation$correlation_stability <- assess_correlation_stability_across_depths(
    simulation_data, numeric_properties, criteria
  )

  return(depth_validation)
}

#' Assess Cholesky Decomposition
#'
#'
#' @param correlation_matrices Correlation matrices and decompositions
#' @param criteria Cholesky validation criteria
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's log level so \code{INFO}-level progress messages print for the duration of this call (default \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @return Cholesky decomposition assessment
#' @keywords internal
diagnose_cholesky <- function(correlation_matrices, criteria = NULL, verbose = getOption("ssurgo.verbose", FALSE)) {

  .old_log_cfg <- set_verbose_logging(verbose)
  on.exit(options(soil_workflow_log_config = .old_log_cfg), add = TRUE)

  if (is.null(criteria)) {
    criteria <- list(
      reconstruction_tolerance = 1e-10,
      condition_number_threshold = 1e12,
      eigenvalue_threshold = 1e-8
    )
  }

  log_message("DEBUG", "Assessing Cholesky decomposition quality", category = "Validation")

  cholesky_validation <- list(
    decomposition_quality = list(),
    numerical_stability = list(),
    reconstruction_accuracy = list()
  )

  # Check different types of correlation matrices with error handling
  matrix_types <- names(correlation_matrices)

  for (matrix_type in matrix_types) {
    matrix_data <- correlation_matrices[[matrix_type]]

    if (is.matrix(matrix_data)) {
      cholesky_validation$decomposition_quality[[matrix_type]] <- tryCatch({
        assess_single_cholesky_decomposition(matrix_data, criteria)
      }, error = function(e) {
        handle_workflow_error(e, paste("Cholesky decomposition assessment for", matrix_type), "warn")
        return(list(assessment_failed = TRUE, error = e$message))
      })
    }
  }

  return(cholesky_validation)
}

# ============================================================================
# 4. GP MODEL VALIDATION
# ============================================================================

#' Validate GP Model Workflow
#'
#'
#' @param gp_models GP models from gp_modeling module
#' @param training_data Original training data
#' @param config GP validation configuration
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's log level so \code{INFO}-level progress messages print for the duration of this call (default \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @return GP model validation results
#' @keywords internal
diagnose_gp_models <- function(gp_models, training_data, config = NULL, verbose = getOption("ssurgo.verbose", FALSE)) {

  .old_log_cfg <- set_verbose_logging(verbose)
  on.exit(options(soil_workflow_log_config = .old_log_cfg), add = TRUE)

  log_message("INFO", "Validating GP model workflow", category = "Validation")

  # Use the shared configuration helpers
  if (is.null(config)) {
    config <- default_config("gp_models")
  }

  # data validation
  if (!is.null(training_data)) {
    training_data_validation <- validate_data_quality(
      training_data,
      required_columns = c("cokey", "hzdept_r"),
      quality_thresholds = config$data_quality_thresholds
    )
  } else {
    training_data_validation <- list(validation_skipped = TRUE)
  }

  validation <- list(
    model_performance = list(),
    prediction_quality = list(),
    depth_trend_realism = list(),
    cross_validation = list(),
    training_data_quality = training_data_validation
  )

  # validation pipeline using the GP depth-modeling functions
  validation_steps <- list(
    "model_performance" = function() diagnose_gp_performance(gp_models, training_data, config$performance_criteria),
    "prediction_quality" = function() diagnose_gp_predictions(gp_models, config$prediction_criteria),
    "depth_trend_realism" = function() diagnose_depth_trend_realism(gp_models, config$realism_criteria),
    "cross_validation" = function() perform_gp_cross_validation(gp_models, training_data, config$cv_criteria)
  )

  for (i in seq_along(validation_steps)) {
    step_name <- names(validation_steps)[i]
    log_message("DEBUG", paste("Validating GP", step_name), category = "Validation")

    track_progress(i, length(validation_steps), "GP model validation", update_frequency = 1)

    validation[[step_name]] <- tryCatch({
      validation_steps[[i]]()
    }, error = function(e) {
      handle_workflow_error(e, paste("GP", step_name, "validation"), "warn")
      return(list(validation_failed = TRUE, error = e$message))
    })
  }

  log_message("INFO", "GP model validation completed", category = "Validation")
  return(validation)
}

#' Validate GP Model Performance
#'
#'
#' @param gp_models GP models
#' @param training_data Training data
#' @param criteria Performance criteria
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's log level so \code{INFO}-level progress messages print for the duration of this call (default \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @return GP model performance assessment
#' @keywords internal
diagnose_gp_performance <- function(gp_models, training_data, criteria = NULL, verbose = getOption("ssurgo.verbose", FALSE)) {

  .old_log_cfg <- set_verbose_logging(verbose)
  on.exit(options(soil_workflow_log_config = .old_log_cfg), add = TRUE)

  if (is.null(criteria)) {
    criteria <- list(
      max_training_rmse = 0.5,
      min_r_squared = 0.6,
      max_condition_number = 1e12
    )
  }

  log_message("DEBUG", "Validating GP model performance", category = "Validation")

  performance_results <- list(
    individual_model_performance = list(),
    overall_performance = list()
  )

  properties <- names(gp_models)[names(gp_models) != "model_summary"]

  for (prop in properties) {
    prop_models <- gp_models[[prop]]

    if (prop_models$type == "stratified_grouped") {
      prop_performance <- list()

      for (group in names(prop_models$models)) {
        group_model <- prop_models$models[[group]]

        if (!is.null(group_model)) {
          prop_performance[[group]] <- tryCatch({
            assess_single_gp_performance(group_model, criteria)
          }, error = function(e) {
            handle_workflow_error(e, paste("GP performance assessment for", prop, group), "warn")
            return(list(assessment_failed = TRUE))
          })
        }
      }

      performance_results$individual_model_performance[[prop]] <- prop_performance
    }
  }

  # Calculate overall performance metrics safe operations
  performance_results$overall_performance <- calculate_overall_gp_performance(
    performance_results$individual_model_performance
  )

  return(performance_results)
}

#' Assess Depth Trend Realism
#'
#'
#' @param gp_models GP models
#' @param criteria Realism criteria
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's log level so \code{INFO}-level progress messages print for the duration of this call (default \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @return Depth trend realism assessment
#' @keywords internal
diagnose_depth_trend_realism <- function(gp_models, criteria = NULL, verbose = getOption("ssurgo.verbose", FALSE)) {

  .old_log_cfg <- set_verbose_logging(verbose)
  on.exit(options(soil_workflow_log_config = .old_log_cfg), add = TRUE)

  if (is.null(criteria)) {
    criteria <- list(
      test_depths = seq(0, 200, by = 5),
      realistic_ranges = get_realistic_property_ranges(),
      monotonicity_tolerance = 0.3,
      gradient_limits = get_realistic_gradients()
    )
  }

  log_message("DEBUG", "Assessing depth trend realism", category = "Validation")

  realism_results <- list(
    trend_predictions = list(),
    realism_assessment = list(),
    constraint_violations = list()
  )

  properties <- names(gp_models)[names(gp_models) != "model_summary"]

  for (prop in properties) {
    prop_models <- gp_models[[prop]]

    if (prop_models$type == "stratified_grouped") {
      prop_realism <- list()

      for (group in names(prop_models$models)) {
        group_model <- prop_models$models[[group]]

        if (!is.null(group_model)) {
          # Use the GP depth-trend predictor
          predictions <- tryCatch({
            predict_gp_depth_trends(group_model, criteria$test_depths)
          }, error = function(e) {
            handle_workflow_error(e, paste("GP prediction for", prop, group), "warn")
            return(NULL)
          })

          if (!is.null(predictions)) {
            # Assess realism validation utilities
            prop_realism[[group]] <- assess_trend_realism(
              predictions, criteria$test_depths, prop, criteria
            )
          }
        }
      }

      realism_results$trend_predictions[[prop]] <- prop_realism
    }
  }

  # Summarize constraint violations patterns
  realism_results$constraint_violations <- summarize_constraint_violations(
    realism_results$trend_predictions
  )

  return(realism_results)
}

#' Validate GP Predictions
#'
#'
#' @param gp_models GP models
#' @param criteria Prediction criteria
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's log level so \code{INFO}-level progress messages print for the duration of this call (default \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @return GP prediction validation
#' @keywords internal
diagnose_gp_predictions <- function(gp_models, criteria = NULL, verbose = getOption("ssurgo.verbose", FALSE)) {

  .old_log_cfg <- set_verbose_logging(verbose)
  on.exit(options(soil_workflow_log_config = .old_log_cfg), add = TRUE)

  if (is.null(criteria)) {
    criteria <- list(
      test_depths = seq(0, 150, by = 10),
      uncertainty_threshold = 0.5,
      smoothness_criteria = 0.1
    )
  }

  log_message("DEBUG", "Validating GP predictions", category = "Validation")

  prediction_validation <- list(
    prediction_smoothness = list(),
    uncertainty_assessment = list(),
    extrapolation_behavior = list()
  )

  properties <- names(gp_models)[names(gp_models) != "model_summary"]

  for (prop in properties) {
    prop_models <- gp_models[[prop]]

    if (prop_models$type == "stratified_grouped") {
      prop_validation <- list()

      for (group in names(prop_models$models)) {
        group_model <- prop_models$models[[group]]

        if (!is.null(group_model)) {
          prop_validation[[group]] <- tryCatch({
            validate_single_gp_predictions(group_model, criteria)
          }, error = function(e) {
            handle_workflow_error(e, paste("GP prediction validation for", prop, group), "warn")
            return(list(validation_failed = TRUE))
          })
        }
      }

      prediction_validation$prediction_smoothness[[prop]] <- prop_validation
    }
  }

  return(prediction_validation)
}

# ============================================================================
# 5. SOIL SCIENCE VALIDATION
# ============================================================================

#' Validate Soil Science Realism
#'
#'
#' @param simulation_data Final simulation data
#' @param original_data Original SSURGO data for comparison
#' @param config Soil science validation configuration
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's log level so \code{INFO}-level progress messages print for the duration of this call (default \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @return Soil science validation results
#' @keywords internal
diagnose_soil_science_realism <- function(simulation_data, original_data = NULL, config = NULL, verbose = getOption("ssurgo.verbose", FALSE)) {

  .old_log_cfg <- set_verbose_logging(verbose)
  on.exit(options(soil_workflow_log_config = .old_log_cfg), add = TRUE)

  log_message("INFO", "Validating soil science realism", category = "Validation")

  # Use the shared configuration helpers
  if (is.null(config)) {
    config <- default_config("soil_science")
  }

  # data validation
  soil_data_validation <- validate_data_quality(
    simulation_data,
    required_columns = c("cokey", "hzdept_r"),
    quality_thresholds = config$data_quality_thresholds
  )

  validation <- list(
    property_constraints = list(),
    horizon_characteristics = list(),
    pedological_relationships = list(),
    depth_trend_realism = list(),
    data_quality = soil_data_validation
  )

  # validation pipeline
  validation_steps <- list(
    "property_constraints" = function() diagnose_property_constraints(simulation_data, config$constraint_criteria),
    "horizon_characteristics" = function() validate_horizon_characteristics(simulation_data, config$horizon_criteria),
    "pedological_relationships" = function() assess_pedological_relationships(simulation_data, config$relationship_criteria),
    "depth_trend_realism" = function() validate_simulation_depth_trends(simulation_data, original_data, config$depth_trend_criteria)
  )

  for (i in seq_along(validation_steps)) {
    step_name <- names(validation_steps)[i]
    log_message("DEBUG", paste("Validating", step_name), category = "Validation")

    track_progress(i, length(validation_steps), "Soil science validation", update_frequency = 1)

    validation[[step_name]] <- tryCatch({
      validation_steps[[i]]()
    }, error = function(e) {
      handle_workflow_error(e, paste("Soil science", step_name, "validation"), "warn")
      return(list(validation_failed = TRUE, error = e$message))
    })
  }

  log_message("INFO", "Soil science validation completed", category = "Validation")
  return(validation)
}

#' Assess Property Constraints
#'
#'
#' @param simulation_data Simulation data
#' @param criteria Property constraint criteria
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's log level so \code{INFO}-level progress messages print for the duration of this call (default \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @return Property constraint assessment
#' @keywords internal
diagnose_property_constraints <- function(simulation_data, criteria = NULL, verbose = getOption("ssurgo.verbose", FALSE)) {

  .old_log_cfg <- set_verbose_logging(verbose)
  on.exit(options(soil_workflow_log_config = .old_log_cfg), add = TRUE)

  if (is.null(criteria)) {
    criteria <- get_default_property_constraints()
  }

  log_message("DEBUG", "Assessing property constraints", category = "Validation")

  constraint_results <- list(
    range_violations = list(),
    cross_property_violations = list(),
    distribution_anomalies = list()
  )

  # property validation
  properties_validation <- validate_properties(names(simulation_data), "laboratory", strict_mode = FALSE)

  # Check range violations range validation
  for (prop in names(criteria$property_ranges)) {
    if (prop %in% names(simulation_data)) {
      range_info <- criteria$property_ranges[[prop]]
      values <- simulation_data[[prop]][!is.na(simulation_data[[prop]])]

      if (length(values) > 0) {
        # Use shared range validation
        range_validation <- validate_numeric_ranges(
          data.frame(temp_prop = values),
          list(temp_prop = range_info),
          action = "warn"
        )

        violations <- values < range_info$min | values > range_info$max

        constraint_results$range_violations[[prop]] <- list(
          n_violations = sum(violations),
          pct_violations = mean(violations) * 100,
          min_value = min(values),
          max_value = max(values),
          expected_range = c(range_info$min, range_info$max),
          module8_validation = range_validation
        )
      }
    }
  }

  # cross-property constraints validation
  constraint_results$cross_property_violations <- assess_cross_property_constraints(
    simulation_data, criteria$cross_property_rules
  )

  # distribution anomaly detection outlier detection
  constraint_results$distribution_anomalies <- detect_distribution_anomalies(
    simulation_data, criteria$distribution_checks
  )

  return(constraint_results)
}

