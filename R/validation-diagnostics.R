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
#' Master validation function enhanced  for comprehensive
#' assessment of the entire soil simulation workflow.
#'
#' @param workflow_results Complete workflow results including all intermediate steps
#' @param original_data Original SSURGO/NRCS data for comparison
#' @param validation_config Configuration for validation parameters (uses Module 0 defaults)
#' @param generate_plots Whether to generate diagnostic plots (default = TRUE)
#' @param output_dir Directory for saving validation outputs
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's log level so \code{INFO}-level progress messages print for the duration of this call (default \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @return Comprehensive validation results
#' @export
validate_complete_workflow <- function(workflow_results,
                                       original_data = NULL,
                                       validation_config = NULL,
                                       generate_plots = TRUE,
                                       output_dir = NULL,
                                       verbose = getOption("ssurgo.verbose", FALSE)) {

  .old_log_cfg <- set_verbose_logging(verbose)
  on.exit(options(soil_workflow_log_config = .old_log_cfg), add = TRUE)

  log_message("INFO", "=== COMPREHENSIVE WORKFLOW VALIDATION ===", category = "Validation")
  start_time <- Sys.time()

  # Use Module 0 configuration management
  if (is.null(validation_config)) {
    validation_config <- get_default_configuration("validation")
  }

  # Initialize validation results structure with Module 0 metadata patterns
  validation_results <- initialize_validation_structure(start_time, validation_config, output_dir)

  # Extract and validate workflow components  utilities
  components <- extract_and_validate_components(workflow_results)

  log_message("INFO", paste("Validating workflow with", length(components), "components"), category = "Validation")

  # Enhanced validation pipeline
  validation_results <- execute_validation_pipeline(
    validation_results, components, original_data, validation_config, generate_plots, output_dir
  )

  # Finalize validation
  validation_results <- finalize_validation_results(validation_results, start_time)

  log_message("INFO", "=== VALIDATION COMPLETE ===", category = "Validation")
  log_message("INFO", paste("Overall quality score:", round(validation_results$overall_assessment$quality_score, 3)), category = "Validation")

  return(validation_results)
}

#' Generate Validation Report
#'
#' Enhanced report generation  I/O utilities and error handling.
#'
#' @param validation_results Results from validate_complete_workflow()
#' @param output_format Format for report: "html", "pdf", or "markdown"
#' @param output_file Output file path
#' @param include_plots Whether to include diagnostic plots
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's log level so \code{INFO}-level progress messages print for the duration of this call (default \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @return Success status
#' @export
generate_validation_report <- function(validation_results,
                                       output_format = "html",
                                       output_file = NULL,
                                       include_plots = TRUE,
                                       verbose = getOption("ssurgo.verbose", FALSE)) {

  .old_log_cfg <- set_verbose_logging(verbose)
  on.exit(options(soil_workflow_log_config = .old_log_cfg), add = TRUE)

  log_message("INFO", "=== GENERATING VALIDATION REPORT ===", category = "Validation")
  log_message("INFO", paste("Output format:", output_format), category = "Validation")

  # Enhanced file naming  utilities
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

  # Create report content with enhanced error handling
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
#' Enhanced quality assessment  validation framework.
#'
#' @param validation_results Complete validation results
#' @param validation_config Validation configuration
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's log level so \code{INFO}-level progress messages print for the duration of this call (default \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @return Overall quality assessment
#' @export
assess_workflow_quality <- function(validation_results, validation_config, verbose = getOption("ssurgo.verbose", FALSE)) {

  .old_log_cfg <- set_verbose_logging(verbose)
  on.exit(options(soil_workflow_log_config = .old_log_cfg), add = TRUE)

  log_message("INFO", "Performing overall workflow quality assessment", category = "Validation")

  # Component quality scores  safe operations
  quality_scores <- calculate_component_quality_scores(validation_results)

  # Calculate weighted overall score  utilities
  weights <- validation_config$quality_weights %||% get_default_quality_weights()
  overall_score <- calculate_weighted_quality_score(quality_scores, weights)

  # Determine quality grade  patterns
  quality_grade <- determine_quality_grade(overall_score)

  # Identify critical issues  validation framework
  critical_issues <- identify_critical_issues(validation_results, validation_config)

  # Create assessment with Module 0 metadata patterns
  assessment <- create_quality_assessment(
    overall_score, quality_grade, quality_scores, critical_issues, validation_config
  )

  log_message("INFO", paste("Quality assessment - Score:", round(overall_score, 3),
                            "Grade:", quality_grade, "Status:", assessment$workflow_status), category = "Validation")

  return(assessment)
}

# ============================================================================
# 2. MONTE CARLO VALIDATION FUNCTIONS (Enhanced)
# ============================================================================

#' Validate Monte Carlo Quality
#'
#' Enhanced Monte Carlo validation  utilities and validation framework.
#'
#' @param monte_carlo_results Results from monte_carlo module
#' @param original_data Original SSURGO data for comparison
#' @param config Monte Carlo validation configuration
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's log level so \code{INFO}-level progress messages print for the duration of this call (default \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @return Monte Carlo validation results
#' @export
validate_monte_carlo_quality <- function(monte_carlo_results,
                                         original_data = NULL,
                                         config = NULL,
                                         verbose = getOption("ssurgo.verbose", FALSE)) {

  .old_log_cfg <- set_verbose_logging(verbose)
  on.exit(options(soil_workflow_log_config = .old_log_cfg), add = TRUE)

  log_message("INFO", "Validating Monte Carlo simulation quality", category = "Validation")

  # Use Module 0 configuration management
  if (is.null(config)) {
    config <- get_default_configuration("monte_carlo")
  }

  # Enhanced data validation
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

  # Enhanced validation pipeline with progress tracking
  validation_steps <- list(
    "convergence" = function() assess_simulation_convergence(simulation_data, config$convergence_criteria),
    "coverage" = function() assess_simulation_coverage(simulation_data, original_data, config$coverage_criteria),
    "distribution" = function() validate_distribution_fidelity(simulation_data, monte_carlo_results$metadata, config$distribution_criteria),
    # generate_simulation_diagnostics() takes 5 args (original_data, simulation_results,
    # properties, correlation_config, config) - a prior 2-arg call here always errored and was
    # silently swallowed by the tryCatch below into a validation_failed placeholder. Sourced from
    # monte_carlo_results, mirroring generate_monte_carlo_realizations()'s own call convention.
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
#' Enhanced coverage assessment  validation utilities.
#'
#' @param simulation_data Simulation results
#' @param original_data Original data for comparison
#' @param criteria Coverage assessment criteria
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's log level so \code{INFO}-level progress messages print for the duration of this call (default \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @return Coverage assessment results
#' @export
assess_simulation_coverage <- function(simulation_data,
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

  # Enhanced property detection
  properties_validation <- validate_properties(names(simulation_data), "laboratory", strict_mode = FALSE)
  numeric_properties <- intersect(
    properties_validation$property_stats$valid_property_names %||% names(simulation_data),
    names(simulation_data)[sapply(simulation_data, is.numeric)]
  )

  # Remove metadata columns
  numeric_properties <- setdiff(numeric_properties, c("simulation_number", "hzdept_r", "hzdepb_r"))

  log_message("DEBUG", paste("Assessing coverage for", length(numeric_properties), "properties"), category = "Validation")

  # Parameter space coverage with enhanced error handling
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

  # Distributional coverage  safe operations
  for (prop in numeric_properties) {
    coverage_results$distributional_coverage[[prop]] <- tryCatch({
      assess_distributional_coverage(simulation_data[[prop]], criteria)
    }, error = function(e) {
      handle_workflow_error(e, paste("Distributional coverage for", prop), "warn")
      return(list(coverage_failed = TRUE))
    })
  }

  # Group representation with Module 0 validation
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
#' Enhanced distribution validation  utilities.
#'
#' @param simulation_data Simulation results
#' @param simulation_metadata Metadata from Monte Carlo generation
#' @param criteria Distribution validation criteria
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's log level so \code{INFO}-level progress messages print for the duration of this call (default \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @return Distribution fidelity assessment
#' @export
validate_distribution_fidelity <- function(simulation_data,
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

  # Enhanced distribution parameter extraction
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

          # Moment comparisons  safe operations
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

          # Outlier assessment  utilities
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
# 3. CORRELATION STRUCTURE VALIDATION (Enhanced with Module 0)
# ============================================================================

#' Validate Correlation Structures
#'
#' Enhanced correlation validation  utilities and Module 6 integration.
#'
#' @param correlation_matrices Correlation matrices from correlation_structure module
#' @param simulation_data Final simulation data
#' @param config Correlation validation configuration
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's log level so \code{INFO}-level progress messages print for the duration of this call (default \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @return Correlation validation results
#' @export
validate_correlation_structures <- function(correlation_matrices,
                                            simulation_data,
                                            config = NULL,
                                            verbose = getOption("ssurgo.verbose", FALSE)) {

  .old_log_cfg <- set_verbose_logging(verbose)
  on.exit(options(soil_workflow_log_config = .old_log_cfg), add = TRUE)

  log_message("INFO", "Validating correlation structures", category = "Validation")

  # Use Module 0 configuration management
  if (is.null(config)) {
    config <- get_default_configuration("correlation")
  }

  # Enhanced data validation
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

  # Enhanced validation pipeline with Module 0 progress tracking
  validation_steps <- list(
    "matrix_quality" = function() validate_correlation_matrix_quality(correlation_matrices, config$matrix_criteria),
    "preservation" = function() validate_correlation_preservation_diagnostics(correlation_matrices, simulation_data, config$preservation_criteria),
    "depth_specific" = function() validate_within_depth_correlations(simulation_data, config$depth_criteria),
    "cholesky" = function() assess_cholesky_decomposition(correlation_matrices, config$cholesky_criteria)
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
#' Enhanced correlation preservation validation using Module 6 utilities.
#'
#' @param original_correlations Original correlation matrices
#' @param simulation_data Final simulation data
#' @param criteria Preservation criteria
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's log level so \code{INFO}-level progress messages print for the duration of this call (default \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @return Correlation preservation assessment
#' @export
validate_correlation_preservation_diagnostics <- function(original_correlations,
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

  # Enhanced property detection
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

  # Overall preservation assessment  safe correlation
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

  # Depth-specific preservation with enhanced error handling
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
#' Enhanced within-depth correlation validation  utilities.
#'
#' @param simulation_data Simulation data
#' @param criteria Depth validation criteria
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's log level so \code{INFO}-level progress messages print for the duration of this call (default \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @return Within-depth correlation validation
#' @export
validate_within_depth_correlations <- function(simulation_data, criteria = NULL, verbose = getOption("ssurgo.verbose", FALSE)) {

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

  # Enhanced property detection
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

  # Assess correlations within depth bins  safe operations
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
#' Enhanced Cholesky validation  error handling.
#'
#' @param correlation_matrices Correlation matrices and decompositions
#' @param criteria Cholesky validation criteria
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's log level so \code{INFO}-level progress messages print for the duration of this call (default \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @return Cholesky decomposition assessment
#' @export
assess_cholesky_decomposition <- function(correlation_matrices, criteria = NULL, verbose = getOption("ssurgo.verbose", FALSE)) {

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

  # Check different types of correlation matrices with enhanced error handling
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
# 4. GP MODEL VALIDATION (Enhanced with Module 5 Integration)
# ============================================================================

#' Validate GP Model Workflow
#'
#' Enhanced GP validation using Module 5 functions and Module 0 utilities.
#'
#' @param gp_models GP models from gp_modeling module
#' @param training_data Original training data
#' @param config GP validation configuration
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's log level so \code{INFO}-level progress messages print for the duration of this call (default \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @return GP model validation results
#' @export
validate_gp_model_workflow <- function(gp_models, training_data, config = NULL, verbose = getOption("ssurgo.verbose", FALSE)) {

  .old_log_cfg <- set_verbose_logging(verbose)
  on.exit(options(soil_workflow_log_config = .old_log_cfg), add = TRUE)

  log_message("INFO", "Validating GP model workflow", category = "Validation")

  # Use Module 0 configuration management
  if (is.null(config)) {
    config <- get_default_configuration("gp_models")
  }

  # Enhanced data validation
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

  # Enhanced validation pipeline using Module 5 functions
  validation_steps <- list(
    "model_performance" = function() validate_gp_model_performance(gp_models, training_data, config$performance_criteria),
    "prediction_quality" = function() validate_gp_predictions(gp_models, config$prediction_criteria),
    "depth_trend_realism" = function() assess_depth_trend_realism(gp_models, config$realism_criteria),
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
#' Enhanced GP performance validation using Module 5 functions.
#'
#' @param gp_models GP models
#' @param training_data Training data
#' @param criteria Performance criteria
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's log level so \code{INFO}-level progress messages print for the duration of this call (default \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @return GP model performance assessment
#' @export
validate_gp_model_performance <- function(gp_models, training_data, criteria = NULL, verbose = getOption("ssurgo.verbose", FALSE)) {

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

  # Calculate overall performance metrics  safe operations
  performance_results$overall_performance <- calculate_overall_gp_performance(
    performance_results$individual_model_performance
  )

  return(performance_results)
}

#' Assess Depth Trend Realism
#'
#' Enhanced depth trend validation using Module 5 prediction functions and Module 0 utilities.
#'
#' @param gp_models GP models
#' @param criteria Realism criteria
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's log level so \code{INFO}-level progress messages print for the duration of this call (default \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @return Depth trend realism assessment
#' @export
assess_depth_trend_realism <- function(gp_models, criteria = NULL, verbose = getOption("ssurgo.verbose", FALSE)) {

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
          # Use Module 5 prediction function
          predictions <- tryCatch({
            predict_gp_depth_trends(group_model, criteria$test_depths)
          }, error = function(e) {
            handle_workflow_error(e, paste("GP prediction for", prop, group), "warn")
            return(NULL)
          })

          if (!is.null(predictions)) {
            # Assess realism  validation utilities
            prop_realism[[group]] <- assess_trend_realism(
              predictions, criteria$test_depths, prop, criteria
            )
          }
        }
      }

      realism_results$trend_predictions[[prop]] <- prop_realism
    }
  }

  # Summarize constraint violations  patterns
  realism_results$constraint_violations <- summarize_constraint_violations(
    realism_results$trend_predictions
  )

  return(realism_results)
}

#' Validate GP Predictions
#'
#' Enhanced GP prediction validation using Module 5 functions.
#'
#' @param gp_models GP models
#' @param criteria Prediction criteria
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's log level so \code{INFO}-level progress messages print for the duration of this call (default \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @return GP prediction validation
#' @export
validate_gp_predictions <- function(gp_models, criteria = NULL, verbose = getOption("ssurgo.verbose", FALSE)) {

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
# 5. SOIL SCIENCE VALIDATION (Enhanced with Module 0)
# ============================================================================

#' Validate Soil Science Realism
#'
#' Enhanced soil science validation  property validation and utilities.
#'
#' @param simulation_data Final simulation data
#' @param original_data Original SSURGO data for comparison
#' @param config Soil science validation configuration
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's log level so \code{INFO}-level progress messages print for the duration of this call (default \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @return Soil science validation results
#' @export
validate_soil_science_realism <- function(simulation_data, original_data = NULL, config = NULL, verbose = getOption("ssurgo.verbose", FALSE)) {

  .old_log_cfg <- set_verbose_logging(verbose)
  on.exit(options(soil_workflow_log_config = .old_log_cfg), add = TRUE)

  log_message("INFO", "Validating soil science realism", category = "Validation")

  # Use Module 0 configuration management
  if (is.null(config)) {
    config <- get_default_configuration("soil_science")
  }

  # Enhanced data validation
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

  # Enhanced validation pipeline with Module 0 progress tracking
  validation_steps <- list(
    "property_constraints" = function() assess_property_constraints(simulation_data, config$constraint_criteria),
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
#' Enhanced property constraint validation  property utilities.
#'
#' @param simulation_data Simulation data
#' @param criteria Property constraint criteria
#' @param verbose Logical; if \code{TRUE}, temporarily raises the package's log level so \code{INFO}-level progress messages print for the duration of this call (default \code{FALSE} - quiet). See \code{set_verbose_logging()}.
#' @return Property constraint assessment
#' @export
assess_property_constraints <- function(simulation_data, criteria = NULL, verbose = getOption("ssurgo.verbose", FALSE)) {

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

  # Enhanced property validation
  properties_validation <- validate_properties(names(simulation_data), "laboratory", strict_mode = FALSE)

  # Check range violations  range validation
  for (prop in names(criteria$property_ranges)) {
    if (prop %in% names(simulation_data)) {
      range_info <- criteria$property_ranges[[prop]]
      values <- simulation_data[[prop]][!is.na(simulation_data[[prop]])]

      if (length(values) > 0) {
        # Use Module 0 range validation
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

  # Enhanced cross-property constraints  validation
  constraint_results$cross_property_violations <- assess_cross_property_constraints(
    simulation_data, criteria$cross_property_rules
  )

  # Enhanced distribution anomaly detection  outlier detection
  constraint_results$distribution_anomalies <- detect_distribution_anomalies(
    simulation_data, criteria$distribution_checks
  )

  return(constraint_results)
}

# ============================================================================
# 6. ENHANCED HELPER FUNCTIONS (Leveraging Module 0)
# ============================================================================

# Enhanced initialization  patterns
initialize_validation_structure <- function(start_time, validation_config, output_dir) {
  list(
    workflow_summary = list(),
    monte_carlo_validation = list(),
    correlation_validation = list(),
    gp_model_validation = list(),
    soil_science_validation = list(),
    performance_metrics = list(),
    overall_assessment = list(),
    diagnostic_plots = list(),
    recommendations = character(0),
    validation_metadata = list(
      timestamp = start_time,
      config = validation_config,
      output_dir = output_dir,
      r_version = R.version.string,
      package_versions = get_package_versions()
    )
  )
}

# Enhanced component extraction with Module 0 validation
extract_and_validate_components <- function(workflow_results) {

  components <- list()

  tryCatch({
    # Enhanced component extraction with validation
    if (is.list(workflow_results)) {

      # Monte Carlo results
      if ("simulation_data" %in% names(workflow_results)) {
        components$monte_carlo_results <- workflow_results
        components$simulation_data <- workflow_results$simulation_data
      }

      # Integrated results with Module 6 structure
      if ("integrated_data" %in% names(workflow_results)) {
        components$final_data <- workflow_results$integrated_data
        components$monte_carlo_results <- list(
          simulation_data = workflow_results$original_simulation_data,
          metadata = workflow_results$original_metadata
        )
      }

      # GP models from Module 5
      if ("gp_models" %in% names(workflow_results)) {
        components$gp_models <- workflow_results$gp_models
      }

      # Correlation matrices
      if ("correlation_matrices" %in% names(workflow_results)) {
        components$correlation_matrices <- workflow_results$correlation_matrices
      }

      # Training data
      if ("training_data" %in% names(workflow_results)) {
        components$training_data <- workflow_results$training_data
      }
    }

    # Validate extracted components
    for (comp_name in names(components)) {
      if (is.data.frame(components[[comp_name]])) {
        validation <- validate_data_quality(
          components[[comp_name]],
          required_columns = character(0),
          quality_thresholds = list()
        )

        if (!validation$overall_quality$validation_passed) {
          log_message("WARN", paste("Component", comp_name, "has data quality issues"), category = "Validation")
        }
      }
    }

  }, error = function(e) {
    handle_workflow_error(e, "Component extraction", "warn")
  })

  return(components)
}

# Enhanced validation pipeline execution
execute_validation_pipeline <- function(validation_results, components, original_data,
                                        validation_config, generate_plots, output_dir) {

  # Validation steps with enhanced error handling
  validation_steps <- list(
    "monte_carlo" = function() {
      if (!is.null(components$monte_carlo_results)) {
        validate_monte_carlo_quality(components$monte_carlo_results, original_data, validation_config$monte_carlo)
      } else {
        log_message("DEBUG", "No Monte Carlo results found", category = "Validation")
        return(list(validation_skipped = TRUE))
      }
    },

    "correlation" = function() {
      if (!is.null(components$correlation_matrices)) {
        validate_correlation_structures(components$correlation_matrices, components$simulation_data, validation_config$correlation)
      } else {
        log_message("DEBUG", "No correlation data found", category = "Validation")
        return(list(validation_skipped = TRUE))
      }
    },

    "gp_models" = function() {
      if (!is.null(components$gp_models)) {
        validate_gp_model_workflow(components$gp_models, components$training_data, validation_config$gp_models)
      } else {
        log_message("DEBUG", "No GP models found", category = "Validation")
        return(list(validation_skipped = TRUE))
      }
    },

    "soil_science" = function() {
      if (!is.null(components$final_data)) {
        validate_soil_science_realism(components$final_data, original_data, validation_config$soil_science)
      } else {
        log_message("DEBUG", "No final simulation data found", category = "Validation")
        return(list(validation_skipped = TRUE))
      }
    }
  )

  # Execute validation steps with progress tracking
  for (i in seq_along(validation_steps)) {
    step_name <- names(validation_steps)[i]
    log_message("INFO", paste("Step", i, ": Validating", step_name), category = "Validation")

    track_progress(i, length(validation_steps) + 3, "Workflow validation", update_frequency = 1)

    validation_results[[paste0(step_name, "_validation")]] <- tryCatch({
      validation_steps[[i]]()
    }, error = function(e) {
      handle_workflow_error(e, paste("Validation step", step_name), "warn")
      return(list(validation_failed = TRUE, error = e$message))
    })
  }

  # Performance metrics
  log_message("INFO", "Step 5: Calculating performance metrics", category = "Validation")
  track_progress(5, 8, "Workflow validation", update_frequency = 1)

  validation_results$performance_metrics <- tryCatch({
    calculate_workflow_performance(components, validation_results)
  }, error = function(e) {
    handle_workflow_error(e, "Performance metrics calculation", "warn")
    return(list(calculation_failed = TRUE))
  })

  # Overall assessment
  log_message("INFO", "Step 6: Performing overall assessment", category = "Validation")
  track_progress(6, 8, "Workflow validation", update_frequency = 1)

  validation_results$overall_assessment <- tryCatch({
    assess_workflow_quality(validation_results, validation_config)
  }, error = function(e) {
    handle_workflow_error(e, "Overall quality assessment", "warn")
    return(list(assessment_failed = TRUE))
  })

  # Generate diagnostic plots
  if (generate_plots) {
    log_message("INFO", "Step 7: Generating diagnostic plots", category = "Validation")
    track_progress(7, 8, "Workflow validation", update_frequency = 1)

    validation_results$diagnostic_plots <- tryCatch({
      generate_comprehensive_diagnostics(components, validation_results, output_dir)
    }, error = function(e) {
      handle_workflow_error(e, "Diagnostic plot generation", "warn")
      return(list(plot_generation_failed = TRUE))
    })
  }

  # Generate recommendations
  log_message("INFO", "Step 8: Generating recommendations", category = "Validation")
  track_progress(8, 8, "Workflow validation", update_frequency = 1)

  validation_results$recommendations <- tryCatch({
    generate_workflow_recommendations(validation_results)
  }, error = function(e) {
    handle_workflow_error(e, "Recommendation generation", "warn")
    return(character(0))
  })

  return(validation_results)
}

# Enhanced finalization with Module 0 metadata
finalize_validation_results <- function(validation_results, start_time) {
  end_time <- Sys.time()
  validation_time <- difftime(end_time, start_time, units = "secs")

  validation_results$validation_metadata$completion_time <- end_time
  validation_results$validation_metadata$validation_duration <- validation_time

  log_message("INFO", paste("Validation time:", round(as.numeric(validation_time), 2), "seconds"), category = "Validation")

  return(validation_results)
}

# Enhanced quality scoring  safe operations
calculate_component_quality_scores <- function(validation_results) {
  quality_scores <- list()

  # Monte Carlo quality  safe operations
  if (!is.null(validation_results$monte_carlo_validation) &&
      !validation_results$monte_carlo_validation$validation_failed %||% FALSE) {
    quality_scores$monte_carlo <- calculate_monte_carlo_score(validation_results$monte_carlo_validation)
  }

  # Correlation preservation quality
  if (!is.null(validation_results$correlation_validation) &&
      !validation_results$correlation_validation$validation_failed %||% FALSE) {
    quality_scores$correlation <- calculate_correlation_score(validation_results$correlation_validation)
  }

  # GP model quality
  if (!is.null(validation_results$gp_model_validation) &&
      !validation_results$gp_model_validation$validation_failed %||% FALSE) {
    quality_scores$gp_models <- calculate_gp_model_score(validation_results$gp_model_validation)
  }

  # Soil science realism quality
  if (!is.null(validation_results$soil_science_validation) &&
      !validation_results$soil_science_validation$validation_failed %||% FALSE) {
    quality_scores$soil_science <- calculate_soil_science_score(validation_results$soil_science_validation)
  }

  return(quality_scores)
}

# Enhanced weighted scoring  safe operations
calculate_weighted_quality_score <- function(quality_scores, weights) {
  if (length(quality_scores) == 0) {
    log_message("WARN", "No quality scores available for weighting", category = "Validation")
    return(0.5)  # Default neutral score
  }

  # Use Module 0 safe operations for weighted mean
  tryCatch({
    score_values <- unlist(quality_scores)
    weight_values <- weights[names(quality_scores)]

    # Handle missing weights
    weight_values[is.na(weight_values)] <- 1.0 / length(quality_scores)

    weighted.mean(score_values, weight_values, na.rm = TRUE)
  }, error = function(e) {
    handle_workflow_error(e, "Weighted quality score calculation", "warn")
    return(mean(unlist(quality_scores), na.rm = TRUE))  # Fallback to simple mean
  })
}

# Enhanced quality grading  patterns
determine_quality_grade <- function(overall_score) {
  if (overall_score >= 0.9) {
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
}

# Enhanced critical issue identification  validation
identify_critical_issues <- function(validation_results, validation_config) {
  critical_issues <- character(0)

  # Check for validation failures
  for (validation_type in c("monte_carlo_validation", "correlation_validation",
                            "gp_model_validation", "soil_science_validation")) {
    if (!is.null(validation_results[[validation_type]]$validation_failed) &&
        validation_results[[validation_type]]$validation_failed) {
      critical_issues <- c(critical_issues, paste("Critical failure in", validation_type))
    }
  }

  # Check for specific quality thresholds
  if (!is.null(validation_results$overall_assessment$quality_score)) {
    if (validation_results$overall_assessment$quality_score < 0.6) {
      critical_issues <- c(critical_issues, "Overall quality score below critical threshold")
    }
  }

  # Check for data quality issues
  for (validation_type in names(validation_results)) {
    if (!is.null(validation_results[[validation_type]]$data_quality) &&
        !validation_results[[validation_type]]$data_quality$overall_quality$validation_passed) {
      critical_issues <- c(critical_issues, paste("Data quality issues in", validation_type))
    }
  }

  return(critical_issues)
}

# Enhanced assessment creation  metadata patterns
create_quality_assessment <- function(overall_score, quality_grade, quality_scores,
                                      critical_issues, validation_config) {
  list(
    quality_score = overall_score,
    quality_grade = quality_grade,
    component_scores = quality_scores,
    critical_issues = critical_issues,
    workflow_status = if (overall_score >= validation_config$minimum_acceptable_score %||% 0.7) {
      "PASSED"
    } else {
      "FAILED"
    },
    confidence_level = calculate_confidence_level(quality_scores),
    assessment_metadata = list(
      timestamp = Sys.time(),
      criteria_used = validation_config$minimum_acceptable_score,
      components_assessed = names(quality_scores)
    )
  )
}

# Enhanced confidence calculation
calculate_confidence_level <- function(quality_scores) {
  if (length(quality_scores) == 0) return(0.5)

  # Calculate confidence based on score consistency and values
  score_values <- unlist(quality_scores)
  mean_score <- mean(score_values, na.rm = TRUE)
  score_variance <- var(score_values, na.rm = TRUE)

  # High mean with low variance = high confidence
  confidence <- mean_score * (1 - min(score_variance, 0.2))

  return(pmax(0, pmin(1, confidence)))
}

## The functions below implement real diagnostics (previously
## placeholder/stub implementations ported as-is from the pre-package
## `modules/` prototype - see [[project_soilsim_package]] stub-completion
## history). `assess_correlation_differences()` just below was already real
## prior to that pass.

#' Assess Simulation Convergence
#'
#' Splits `simulation_data`'s realizations (`simulation_number`) into
#' sequential batches and checks whether each numeric property's running
#' mean has stabilized (relative change between the last two batches below
#' `criteria$tolerance`).
#'
#' @param simulation_data Monte Carlo simulation data (must have
#'   `simulation_number` and numeric property columns).
#' @param criteria List, optionally with `tolerance` (default 0.05) and
#'   `n_batches` (default 5).
#' @return List with `converged`, `convergence_metric` (max relative change
#'   across assessed properties), `assessment_method`.
assess_simulation_convergence <- function(simulation_data, criteria) {
  tolerance <- criteria$tolerance %||% 0.05
  n_batches <- criteria$n_batches %||% 5

  if (is.null(simulation_data) || !("simulation_number" %in% names(simulation_data))) {
    return(list(converged = NA, convergence_metric = NA_real_, assessment_method = "batch_mean_stabilization"))
  }

  numeric_cols <- names(simulation_data)[vapply(simulation_data, is.numeric, logical(1))]
  numeric_cols <- setdiff(numeric_cols, c("simulation_number", "hzdept_r", "hzdepb_r"))

  sim_numbers <- sort(unique(simulation_data$simulation_number))
  n_batches <- max(2, min(n_batches, length(sim_numbers)))
  batch_id <- cut(seq_along(sim_numbers), breaks = n_batches, labels = FALSE)
  batch_of <- stats::setNames(batch_id, sim_numbers)

  relative_changes <- c()
  for (col in numeric_cols) {
    batch_means <- vapply(seq_len(n_batches), function(b) {
      keep <- simulation_data$simulation_number %in% sim_numbers[batch_of == b]
      mean(simulation_data[[col]][keep], na.rm = TRUE)
    }, numeric(1))

    if (length(batch_means) >= 2 && all(is.finite(batch_means[(n_batches - 1):n_batches]))) {
      last_two <- batch_means[(n_batches - 1):n_batches]
      denom <- max(abs(last_two[1]), 1e-8)
      relative_changes <- c(relative_changes, abs(diff(last_two)) / denom)
    }
  }

  convergence_metric <- if (length(relative_changes) > 0) max(relative_changes) else NA_real_

  list(
    converged = if (is.na(convergence_metric)) NA else convergence_metric <= tolerance,
    convergence_metric = convergence_metric,
    assessment_method = "batch_mean_stabilization"
  )
}

#' Assess Property Coverage
#'
#' Compares the range/quantile coverage of `sim_values` against `orig_values`:
#' what fraction of `orig_values`' range falls within `sim_values`' range, and
#' whether `sim_values` extrapolates meaningfully beyond it.
#'
#' @param sim_values Simulated values.
#' @param orig_values Original/reference values.
#' @param criteria List, optionally with `max_extrapolation_factor` (default 1.2).
#' @return List with `coverage_percentage`, `extrapolation_detected`, `coverage_quality`.
assess_property_coverage <- function(sim_values, orig_values, criteria) {
  sim_values <- sim_values[is.finite(sim_values)]
  orig_values <- orig_values[is.finite(orig_values)]
  extrapolation_factor <- criteria$max_extrapolation_factor %||% 1.2

  if (length(sim_values) == 0 || length(orig_values) == 0) {
    return(list(coverage_percentage = NA_real_, extrapolation_detected = NA, coverage_quality = "unknown"))
  }

  orig_range <- range(orig_values)
  sim_range <- range(sim_values)

  overlap_lo <- max(orig_range[1], sim_range[1])
  overlap_hi <- min(orig_range[2], sim_range[2])
  overlap <- max(0, overlap_hi - overlap_lo)
  orig_span <- diff(orig_range)

  coverage_percentage <- if (orig_span > 0) (overlap / orig_span) * 100 else 100

  orig_center <- mean(orig_range)
  allowed_half_span <- (orig_span / 2) * extrapolation_factor
  extrapolation_detected <- (sim_range[1] < orig_center - allowed_half_span) ||
    (sim_range[2] > orig_center + allowed_half_span)

  list(
    coverage_percentage = coverage_percentage,
    extrapolation_detected = extrapolation_detected,
    coverage_quality = if (coverage_percentage >= 90 && !extrapolation_detected) "good" else "poor"
  )
}

#' Assess Distributional Coverage
#'
#' Assesses `values`' own distributional spread: outlier percentage (IQR
#' method) and whether the interquartile range is non-degenerate.
#'
#' @param values Values to assess.
#' @param criteria List, optionally with `max_outlier_percentage` (default 10).
#' @return List with `distribution_coverage`, `outlier_percentage`, `coverage_adequate`.
assess_distributional_coverage <- function(values, criteria) {
  values <- values[is.finite(values)]
  max_outlier_pct <- criteria$max_outlier_percentage %||% 10

  if (length(values) < 4) {
    return(list(distribution_coverage = NA_real_, outlier_percentage = NA_real_, coverage_adequate = NA))
  }

  outliers <- tryCatch(detect_outliers(values, method = "iqr", return_indices = FALSE),
                       error = function(e) rep(FALSE, length(values)))
  outlier_percentage <- mean(outliers, na.rm = TRUE) * 100

  q <- stats::quantile(values, probs = c(0.05, 0.95), na.rm = TRUE, names = FALSE)
  full_range <- diff(range(values))
  distribution_coverage <- if (full_range > 0) diff(q) / full_range else 0

  list(
    distribution_coverage = distribution_coverage,
    outlier_percentage = outlier_percentage,
    coverage_adequate = outlier_percentage <= max_outlier_pct
  )
}

#' Assess Group Representation
#'
#' Assesses whether groups (`cokey`) in `simulation_data` are adequately
#' represented (each has at least `criteria$min_samples_per_group` rows).
#'
#' @param simulation_data Simulation data with a `cokey` column.
#' @param criteria List, optionally with `min_samples_per_group` (default 10).
#' @return List with `groups_well_represented` (fraction), `underrepresented_groups`,
#'   `representation_adequate`.
assess_group_representation <- function(simulation_data, criteria) {
  min_samples <- criteria$min_samples_per_group %||% 10

  if (is.null(simulation_data) || !("cokey" %in% names(simulation_data))) {
    return(list(groups_well_represented = NA_real_, underrepresented_groups = character(0),
                representation_adequate = NA))
  }

  group_counts <- table(simulation_data$cokey)
  adequate <- group_counts >= min_samples
  underrepresented <- names(group_counts)[!adequate]

  list(
    groups_well_represented = if (length(group_counts) > 0) mean(adequate) else NA_real_,
    underrepresented_groups = underrepresented,
    representation_adequate = length(underrepresented) == 0
  )
}

#' Test Distribution Fidelity
#'
#' Compares `simulated_values` against a reference sample drawn from the
#' fitted distribution (`param_info$family`/`param_info$fit`, via
#' [quantile_from_fit()]) using a two-sample Kolmogorov-Smirnov test.
#'
#' @param simulated_values Simulated values.
#' @param param_info Fitted distribution parameters, with `family` and `fit`.
#' @param criteria List, optionally with `ks_test_alpha` (default 0.05).
#' @return List with `ks_test_p_value`, `distribution_match`, `test_method`.
test_distribution_fidelity <- function(simulated_values, param_info, criteria) {
  simulated_values <- simulated_values[is.finite(simulated_values)]
  alpha <- criteria$ks_test_alpha %||% 0.05

  if (length(simulated_values) < 2 || is.null(param_info$family) || is.null(param_info$fit)) {
    return(list(ks_test_p_value = NA_real_, distribution_match = NA, test_method = "ks_two_sample"))
  }

  reference_probs <- stats::ppoints(max(length(simulated_values), 200))
  reference_sample <- quantile_from_fit(reference_probs, param_info$family, param_info$fit)
  reference_sample <- reference_sample[is.finite(reference_sample)]

  if (length(reference_sample) < 2) {
    return(list(ks_test_p_value = NA_real_, distribution_match = NA, test_method = "ks_two_sample"))
  }

  ks_result <- tryCatch(stats::ks.test(simulated_values, reference_sample), error = function(e) NULL)
  p_value <- if (!is.null(ks_result)) ks_result$p.value else NA_real_

  list(
    ks_test_p_value = p_value,
    distribution_match = if (is.na(p_value)) NA else p_value >= alpha,
    test_method = "ks_two_sample"
  )
}

#' Compare Distribution Moments
#'
#' Compares `simulated_values`' mean/variance against a reference sample
#' drawn from the fitted distribution.
#'
#' @param simulated_values Simulated values.
#' @param param_info Fitted distribution parameters, with `family` and `fit`.
#' @param criteria List, optionally with `moment_tolerance` (default 0.1).
#' @return List with `mean_difference`, `variance_difference` (relative),
#'   `moment_match_quality`.
compare_distribution_moments <- function(simulated_values, param_info, criteria) {
  simulated_values <- simulated_values[is.finite(simulated_values)]
  tolerance <- criteria$moment_tolerance %||% 0.1

  if (length(simulated_values) < 2 || is.null(param_info$family) || is.null(param_info$fit)) {
    return(list(mean_difference = NA_real_, variance_difference = NA_real_, moment_match_quality = "unknown"))
  }

  reference_probs <- stats::ppoints(max(length(simulated_values), 500))
  reference_sample <- quantile_from_fit(reference_probs, param_info$family, param_info$fit)
  reference_sample <- reference_sample[is.finite(reference_sample)]

  if (length(reference_sample) < 2) {
    return(list(mean_difference = NA_real_, variance_difference = NA_real_, moment_match_quality = "unknown"))
  }

  mean_diff <- abs(mean(simulated_values) - mean(reference_sample)) / max(abs(mean(reference_sample)), 1e-8)
  var_diff <- abs(stats::var(simulated_values) - stats::var(reference_sample)) / max(abs(stats::var(reference_sample)), 1e-8)

  list(
    mean_difference = mean_diff,
    variance_difference = var_diff,
    moment_match_quality = if (mean_diff <= tolerance && var_diff <= tolerance) "good" else "poor"
  )
}

#' Compare Distribution Quantiles
#'
#' Compares `simulated_values`' 10th/50th/90th percentiles against the
#' fitted distribution's theoretical quantiles at those same probabilities.
#'
#' @param simulated_values Simulated values.
#' @param param_info Fitted distribution parameters, with `family` and `fit`.
#' @param criteria List, optionally with `quantile_tolerance` (default 0.05).
#' @return List with `quantile_differences` (relative, length 3),
#'   `quantile_match_quality`.
compare_distribution_quantiles <- function(simulated_values, param_info, criteria) {
  simulated_values <- simulated_values[is.finite(simulated_values)]
  tolerance <- criteria$quantile_tolerance %||% 0.05
  probs <- c(0.1, 0.5, 0.9)

  if (length(simulated_values) < 2 || is.null(param_info$family) || is.null(param_info$fit)) {
    return(list(quantile_differences = rep(NA_real_, length(probs)), quantile_match_quality = "unknown"))
  }

  theoretical <- quantile_from_fit(probs, param_info$family, param_info$fit)
  empirical <- stats::quantile(simulated_values, probs = probs, na.rm = TRUE, names = FALSE)

  quantile_differences <- abs(empirical - theoretical) / pmax(abs(theoretical), 1e-8)

  list(
    quantile_differences = quantile_differences,
    quantile_match_quality = if (all(is.finite(quantile_differences)) && all(quantile_differences <= tolerance)) "good" else "poor"
  )
}

#' Validate Correlation Matrix Quality
#'
#' Inspects every matrix in `correlation_matrices` (a named list) for
#' positive-definiteness (via eigenvalues) and condition number, reusing
#' [assess_correlation_matrix_properties()] per matrix and aggregating
#' across all of them.
#'
#' @param correlation_matrices Named list of correlation matrices.
#' @param criteria List, optionally with `condition_number_threshold`
#'   (default 1e12).
#' @return List with `positive_definite` (TRUE only if all matrices are),
#'   `condition_number` (max across matrices), `matrix_quality`.
validate_correlation_matrix_quality <- function(correlation_matrices, criteria) {
  threshold <- criteria$condition_number_threshold %||% 1e12

  matrices <- Filter(is.matrix, correlation_matrices)
  if (length(matrices) == 0) {
    return(list(positive_definite = NA, condition_number = NA_real_, matrix_quality = "unknown"))
  }

  per_matrix <- lapply(matrices, assess_correlation_matrix_properties)
  all_pd <- all(vapply(per_matrix, function(x) isTRUE(x$is_positive_definite), logical(1)))
  max_condition <- suppressWarnings(max(vapply(per_matrix, function(x) x$condition_number, numeric(1)), na.rm = TRUE))

  list(
    positive_definite = all_pd,
    condition_number = max_condition,
    matrix_quality = if (all_pd && is.finite(max_condition) && max_condition <= threshold) "good" else "poor"
  )
}

#' Assess Correlation Differences
#'
#' Unlike most of its neighbors in this section, this one performs real
#' computation: element-wise absolute difference between two correlation
#' matrices, summarized as max/mean/RMSE.
#'
#' @param orig_cor Original correlation matrix.
#' @param sim_cor Simulated correlation matrix (same dimensions as `orig_cor`).
#' @param criteria Unused.
#' @return List with `max_difference`, `mean_difference`, `correlation_rmse`.
assess_correlation_differences <- function(orig_cor, sim_cor, criteria) {
  cor_diff <- abs(orig_cor - sim_cor)
  list(
    max_difference = max(cor_diff, na.rm = TRUE),
    mean_difference = mean(cor_diff, na.rm = TRUE),
    correlation_rmse = sqrt(mean(cor_diff^2, na.rm = TRUE)),
    assessment_quality = "enhanced"
  )
}

#' Assess Depth-Specific Correlation Preservation
#'
#' Compares a depth-specific correlation matrix (`depth_cor`) against the
#' overall/original correlation structure, reusing
#' [assess_correlation_differences()] on the common properties.
#'
#' @param depth_cor Correlation matrix at one depth.
#' @param original_correlations List with a `global_correlation_matrix` entry.
#' @param depth Depth value (used only for logging/context).
#' @param criteria Assessment criteria (passed through to
#'   [assess_correlation_differences()]).
#' @return List with `preservation_quality` (1 - mean difference, clamped to
#'   `[0, 1]`), `correlation_differences` (mean absolute difference),
#'   `depth_specific_assessment`.
assess_depth_correlation_preservation <- function(depth_cor, original_correlations, depth, criteria) {
  original_cor <- original_correlations$global_correlation_matrix

  if (is.null(original_cor) || is.null(depth_cor) || !is.matrix(depth_cor)) {
    return(list(preservation_quality = NA_real_, correlation_differences = NA_real_,
                depth_specific_assessment = FALSE))
  }

  common_props <- intersect(rownames(original_cor), rownames(depth_cor))
  if (length(common_props) < 2) {
    return(list(preservation_quality = NA_real_, correlation_differences = NA_real_,
                depth_specific_assessment = FALSE))
  }

  diffs <- assess_correlation_differences(
    original_cor[common_props, common_props], depth_cor[common_props, common_props], criteria
  )

  list(
    preservation_quality = max(0, 1 - diffs$mean_difference),
    correlation_differences = diffs$mean_difference,
    depth_specific_assessment = TRUE
  )
}

#' Assess Correlation Matrix Properties
#'
#' Checks a correlation matrix's positive-definiteness (via eigenvalues) and
#' condition number.
#'
#' @param cor_matrix Correlation matrix.
#' @return List with `is_positive_definite`, `condition_number`,
#'   `eigenvalue_range`, `matrix_stability`.
assess_correlation_matrix_properties <- function(cor_matrix) {
  if (is.null(cor_matrix) || !is.matrix(cor_matrix) || nrow(cor_matrix) != ncol(cor_matrix)) {
    return(list(is_positive_definite = NA, condition_number = NA_real_,
                eigenvalue_range = c(NA_real_, NA_real_), matrix_stability = "unknown"))
  }

  eigenvalues <- tryCatch(eigen(cor_matrix, only.values = TRUE, symmetric = TRUE)$values,
                          error = function(e) NA_real_)

  if (all(is.na(eigenvalues))) {
    return(list(is_positive_definite = NA, condition_number = NA_real_,
                eigenvalue_range = c(NA_real_, NA_real_), matrix_stability = "unknown"))
  }

  is_pd <- all(eigenvalues > -1e-8)
  positive_eigenvalues <- eigenvalues[eigenvalues > 1e-12]
  condition_number <- if (length(positive_eigenvalues) > 0) {
    max(eigenvalues) / min(positive_eigenvalues)
  } else {
    Inf
  }

  list(
    is_positive_definite = is_pd,
    condition_number = condition_number,
    eigenvalue_range = range(eigenvalues),
    matrix_stability = if (is_pd && is.finite(condition_number) && condition_number < 1e12) "good" else "poor"
  )
}

#' Assess Correlation Stability Across Depths
#'
#' Computes a correlation matrix for `numeric_properties` within each depth
#' bin of `simulation_data`, and measures how stable pairwise correlations
#' are across bins (1 - mean standard deviation of each pairwise
#' correlation across bins).
#'
#' @param simulation_data Simulation data with an `hzdept_r` column.
#' @param numeric_properties Character vector of property columns to correlate.
#' @param criteria List, optionally with `depth_bins` (default 5 quantile bins).
#' @return List with `stability_score`, `unstable_depths`, `stability_assessment`.
assess_correlation_stability_across_depths <- function(simulation_data, numeric_properties, criteria) {
  if (is.null(simulation_data) || !("hzdept_r" %in% names(simulation_data)) ||
      length(numeric_properties) < 2) {
    return(list(stability_score = NA_real_, unstable_depths = character(0), stability_assessment = "unknown"))
  }

  depth_bins <- criteria$depth_bins %||% {
    breaks <- unique(stats::quantile(simulation_data$hzdept_r, probs = seq(0, 1, length.out = 6), na.rm = TRUE))
    if (length(breaks) < 3) NULL else breaks
  }

  if (is.null(depth_bins)) {
    return(list(stability_score = NA_real_, unstable_depths = character(0), stability_assessment = "unknown"))
  }

  bin <- cut(simulation_data$hzdept_r, breaks = depth_bins, include.lowest = TRUE)
  bin_cors <- list()

  for (b in levels(bin)) {
    bin_data <- simulation_data[!is.na(bin) & bin == b, numeric_properties, drop = FALSE]
    if (nrow(bin_data) >= 5) {
      bin_cor <- tryCatch(stats::cor(bin_data, use = "pairwise.complete.obs"), error = function(e) NULL)
      if (!is.null(bin_cor)) bin_cors[[b]] <- bin_cor
    }
  }

  if (length(bin_cors) < 2) {
    return(list(stability_score = NA_real_, unstable_depths = character(0), stability_assessment = "unknown"))
  }

  pair_idx <- which(upper.tri(bin_cors[[1]]))
  pair_values <- vapply(bin_cors, function(m) m[pair_idx], numeric(length(pair_idx)))
  pair_sd <- apply(pair_values, 1, stats::sd, na.rm = TRUE)

  stability_score <- max(0, 1 - mean(pair_sd, na.rm = TRUE))
  unstable_bins <- names(bin_cors)[apply(pair_values, 2, function(col) any(abs(col - rowMeans(pair_values, na.rm = TRUE)) > 0.3))]

  list(
    stability_score = stability_score,
    unstable_depths = unstable_bins,
    stability_assessment = if (stability_score >= 0.7) "good" else "poor"
  )
}

#' Assess a Single Cholesky Decomposition
#'
#' Actually attempts a Cholesky decomposition of `matrix_data` and measures
#' reconstruction error (`||L'L - matrix_data||` in Frobenius norm).
#'
#' @param matrix_data Matrix to decompose.
#' @param criteria List, optionally with `reconstruction_tolerance` (default 1e-10).
#' @return List with `decomposition_successful`, `reconstruction_error`,
#'   `numerical_stability`.
assess_single_cholesky_decomposition <- function(matrix_data, criteria) {
  tolerance <- criteria$reconstruction_tolerance %||% 1e-10

  chol_result <- tryCatch(chol(matrix_data), error = function(e) NULL)

  if (is.null(chol_result)) {
    return(list(decomposition_successful = FALSE, reconstruction_error = NA_real_,
                numerical_stability = "poor"))
  }

  reconstruction_error <- norm(t(chol_result) %*% chol_result - matrix_data, type = "F")

  list(
    decomposition_successful = TRUE,
    reconstruction_error = reconstruction_error,
    numerical_stability = if (reconstruction_error <= tolerance) "excellent" else if (reconstruction_error <= 1e-6) "good" else "poor"
  )
}

# Enhanced GP assessment functions

#' Assess a Single GP Model's Performance
#'
#' Delegates to the already-real `calculate_model_diagnostics()`
#' (`gp-modeling.R`) for training RMSE, and derives an R-squared from that
#' RMSE against the training data's own variance.
#'
#' @param group_model Fitted GP model (as produced by
#'   `fit_individual_gp_model()`: `gp_model`, `training_data`, `property`).
#' @param criteria List, optionally with `max_training_rmse` (default `Inf`)
#'   and `min_r_squared` (default 0).
#' @return List with `training_rmse`, `r_squared`, `model_quality`,
#'   `performance_assessment`.
assess_single_gp_performance <- function(group_model, criteria) {
  max_rmse <- criteria$max_training_rmse %||% Inf
  min_r2 <- criteria$min_r_squared %||% 0

  if (is.null(group_model) || is.null(group_model$gp_model) || is.null(group_model$training_data)) {
    return(list(training_rmse = NA_real_, r_squared = NA_real_, model_quality = "unknown",
                performance_assessment = "insufficient_data"))
  }

  diagnostics <- calculate_model_diagnostics(group_model$gp_model, group_model$training_data,
                                             group_model$property %||% "unknown")

  training_rmse <- diagnostics$training_rmse
  total_var <- stats::var(group_model$training_data$mean_value, na.rm = TRUE)

  r_squared <- if (is.finite(training_rmse) && is.finite(total_var) && total_var > 0) {
    max(0, 1 - (training_rmse^2 / total_var))
  } else {
    NA_real_
  }

  quality <- if (is.finite(training_rmse) && is.finite(r_squared) &&
                 training_rmse <= max_rmse && r_squared >= min_r2) "good" else "poor"

  list(
    training_rmse = training_rmse,
    r_squared = r_squared,
    model_quality = quality,
    performance_assessment = "computed_from_training_residuals"
  )
}

#' Calculate Overall GP Performance
#'
#' Aggregates the per-property, per-group results of
#' [assess_single_gp_performance()] into overall mean RMSE/R-squared.
#'
#' @param individual_performance Nested list: property -> group -> performance
#'   result (as returned by [assess_single_gp_performance()]).
#' @return List with `mean_r_squared`, `mean_rmse`, `overall_quality`,
#'   `performance_summary`.
calculate_overall_gp_performance <- function(individual_performance) {
  all_r2 <- c()
  all_rmse <- c()

  for (prop in names(individual_performance)) {
    for (group in names(individual_performance[[prop]])) {
      perf <- individual_performance[[prop]][[group]]
      if (!is.null(perf$r_squared) && is.finite(perf$r_squared)) all_r2 <- c(all_r2, perf$r_squared)
      if (!is.null(perf$training_rmse) && is.finite(perf$training_rmse)) all_rmse <- c(all_rmse, perf$training_rmse)
    }
  }

  mean_r2 <- if (length(all_r2) > 0) mean(all_r2) else NA_real_
  mean_rmse <- if (length(all_rmse) > 0) mean(all_rmse) else NA_real_

  list(
    mean_r_squared = mean_r2,
    mean_rmse = mean_rmse,
    overall_quality = if (is.na(mean_r2)) "unknown" else if (mean_r2 >= 0.6) "good" else "poor",
    performance_summary = paste(length(all_r2), "GP model(s) assessed")
  )
}

#' Assess Depth-Trend Realism
#'
#' Reuses the already-real `assess_trend_monotonicity()` and
#' `assess_realistic_values()` (`gp-modeling.R`) to judge whether a
#' predicted depth trend is realistic, and counts constraint violations
#' against `criteria$realistic_ranges` when supplied.
#'
#' @param predictions GP predictions at `test_depths`.
#' @param test_depths Depths predictions were made at.
#' @param prop Property name (matched against `criteria$realistic_ranges`
#'   and `assess_realistic_values()`'s built-in ranges).
#' @param criteria List, optionally with `realistic_ranges` (named list of
#'   `list(min=, max=)`).
#' @return List with `realistic_trend`, `constraint_violations`,
#'   `trend_quality`, `realism_assessment`.
assess_trend_realism <- function(predictions, test_depths, prop, criteria) {
  monotonic <- assess_trend_monotonicity(predictions)
  values_realistic <- assess_realistic_values(predictions, prop)

  realistic_ranges <- criteria$realistic_ranges %||% list()
  violations <- if (prop %in% names(realistic_ranges)) {
    range_vals <- realistic_ranges[[prop]]
    sum(predictions < range_vals$min | predictions > range_vals$max, na.rm = TRUE)
  } else {
    sum(!is.finite(predictions))
  }

  realistic_trend <- isTRUE(values_realistic) && (is.na(monotonic) || isTRUE(monotonic))

  list(
    realistic_trend = realistic_trend,
    constraint_violations = violations,
    trend_quality = if (realistic_trend) "good" else "poor",
    realism_assessment = "computed_from_monotonicity_and_range_checks"
  )
}

#' Summarize Constraint Violations
#'
#' Aggregates the nested property -> group -> [assess_trend_realism()]
#' results into an overall violation count and compliance rate.
#'
#' @param trend_predictions Nested list: property -> group -> trend
#'   assessment (as returned by [assess_trend_realism()]).
#' @return List with `total_violations`, `violation_types` (properties with
#'   at least one non-realistic trend), `overall_compliance`.
summarize_constraint_violations <- function(trend_predictions) {
  total_violations <- 0
  violation_types <- character(0)
  n_assessments <- 0
  n_realistic <- 0

  for (prop in names(trend_predictions)) {
    for (group in names(trend_predictions[[prop]])) {
      result <- trend_predictions[[prop]][[group]]
      if (!is.null(result$constraint_violations)) {
        total_violations <- total_violations + result$constraint_violations
        n_assessments <- n_assessments + 1
        if (isTRUE(result$realistic_trend)) {
          n_realistic <- n_realistic + 1
        } else {
          violation_types <- union(violation_types, prop)
        }
      }
    }
  }

  list(
    total_violations = total_violations,
    violation_types = violation_types,
    overall_compliance = if (n_assessments > 0) n_realistic / n_assessments else NA_real_
  )
}

#' Validate a Single GP Model's Predictions
#'
#' Predicts a depth trend from `group_model` and measures its smoothness as
#' one minus the normalized mean absolute second difference of the
#' predicted curve (a roughness measure).
#'
#' @param group_model Fitted GP model.
#' @param criteria List, optionally with `test_depths` (default
#'   `seq(0, 150, by = 10)`) and `smoothness_criteria` (default 0.1).
#' @return List with `prediction_smoothness`, `uncertainty_reasonable`,
#'   `prediction_quality`.
validate_single_gp_predictions <- function(group_model, criteria) {
  test_depths <- criteria$test_depths %||% seq(0, 150, by = 10)
  smoothness_threshold <- criteria$smoothness_criteria %||% 0.1

  predictions <- tryCatch(predict_gp_depth_trends(group_model, test_depths),
                          error = function(e) rep(NA_real_, length(test_depths)))
  predictions <- predictions[is.finite(predictions)]

  if (length(predictions) < 3) {
    return(list(prediction_smoothness = NA_real_, uncertainty_reasonable = NA,
                prediction_quality = "insufficient_data"))
  }

  roughness <- mean(abs(diff(predictions, differences = 2)))
  value_range <- max(predictions) - min(predictions)
  normalized_roughness <- if (value_range > 0) roughness / value_range else 0
  smoothness <- max(0, 1 - normalized_roughness)

  list(
    prediction_smoothness = smoothness,
    uncertainty_reasonable = smoothness >= (1 - smoothness_threshold),
    prediction_quality = if (smoothness >= (1 - smoothness_threshold)) "good" else "poor"
  )
}

# Enhanced soil science functions
get_default_property_constraints <- function() {
  # NOTE: this previously called a nonexistent get_default_property_constraints()
  # (a dangling reference - the function was never implemented anywhere in
  # modules/), which meant assess_property_constraints() always errored with
  # "could not find function" whenever called with criteria = NULL. Fixed by
  # building constraints$property_ranges from the real, already-implemented
  # get_realistic_property_ranges(), which returns exactly the
  # list(propname = list(min=, max=), ...) shape this function's caller expects.
  constraints <- list(property_ranges = get_realistic_property_ranges())

  # Add enhanced constraints
  constraints$enhanced_checks <- list(
    use_module8_validation = TRUE,
    apply_soil_specific_rules = TRUE
  )

  return(constraints)
}

#' Assess Cross-Property Constraints
#'
#' Checks the texture sum-to-100 constraint (when `sandtotal`/`claytotal`/
#' `silttotal` are present) and any explicit `cross_property_rules` (in the
#' same shape as [create_validation_config()]'s `relationship_rules` -
#' `list(properties=, type=, expected_sum=, tolerance=)`).
#'
#' @param simulation_data Simulation data.
#' @param cross_property_rules List of relationship rules, or `NULL`.
#' @return List with `texture_sum_violations`, `relationship_violations`,
#'   `overall_compliance`.
assess_cross_property_constraints <- function(simulation_data, cross_property_rules) {
  texture_cols <- intersect(c("sandtotal", "claytotal", "silttotal"), names(simulation_data))
  texture_sum_violations <- 0

  if (length(texture_cols) == 3) {
    texture_sum <- rowSums(simulation_data[texture_cols], na.rm = TRUE)
    texture_sum_violations <- sum(abs(texture_sum - 100) > 5, na.rm = TRUE)
  }

  relationship_violations <- 0
  for (rule in cross_property_rules) {
    if (identical(rule$type, "sum") && all(rule$properties %in% names(simulation_data))) {
      rule_sum <- rowSums(simulation_data[rule$properties], na.rm = TRUE)
      expected <- rule$expected_sum %||% 100
      tol <- rule$tolerance %||% 0.1
      relationship_violations <- relationship_violations +
        sum(abs(rule_sum - expected) > (abs(expected) * tol), na.rm = TRUE)
    }
  }

  n_rows <- max(nrow(simulation_data), 1)
  overall_compliance <- 1 - min(1, (texture_sum_violations + relationship_violations) / n_rows)

  list(
    texture_sum_violations = texture_sum_violations,
    relationship_violations = relationship_violations,
    overall_compliance = overall_compliance
  )
}

#' Detect Distribution Anomalies
#'
#' Runs IQR-based outlier detection (reusing the already-real
#' [detect_outliers()]) across `simulation_data`'s numeric property columns.
#'
#' @param simulation_data Simulation data.
#' @param distribution_checks Unused (no per-check configuration is
#'   currently defined upstream - reserved for future per-property checks).
#' @return List with `anomalies_detected` (count), `anomaly_types`,
#'   `anomaly_severity`.
detect_distribution_anomalies <- function(simulation_data, distribution_checks) {
  numeric_cols <- names(simulation_data)[vapply(simulation_data, is.numeric, logical(1))]
  numeric_cols <- setdiff(numeric_cols, c("simulation_number", "hzdept_r", "hzdepb_r"))

  anomalies_detected <- 0
  anomaly_types <- character(0)

  for (col in numeric_cols) {
    outliers <- tryCatch(detect_outliers(simulation_data[[col]], method = "iqr", return_indices = FALSE),
                         error = function(e) rep(FALSE, nrow(simulation_data)))
    n_out <- sum(outliers, na.rm = TRUE)
    if (n_out > 0) {
      anomalies_detected <- anomalies_detected + n_out
      anomaly_types <- union(anomaly_types, "outliers")
    }
  }

  n_rows <- max(nrow(simulation_data), 1)
  severity <- if (anomalies_detected == 0) "none" else if (anomalies_detected / n_rows < 0.05) "minor" else "major"

  list(
    anomalies_detected = anomalies_detected,
    anomaly_types = anomaly_types,
    anomaly_severity = severity
  )
}

#' Validate Horizon Characteristics
#'
#' Checks whether surface (`hzdept_r <= 30`) and subsurface horizons fall
#' within realistic property ranges (`criteria$realistic_ranges`, defaulting
#' to `get_realistic_property_ranges()`), and measures the smoothness of the
#' depth transition for the first available property.
#'
#' @param simulation_data Simulation data.
#' @param criteria List, optionally with `realistic_ranges`.
#' @return List with `surface_horizon_quality`, `subsurface_quality`,
#'   `transition_quality` (each a fraction/score in `[0, 1]`, or `NA` when
#'   not computable).
validate_horizon_characteristics <- function(simulation_data, criteria) {
  ranges <- criteria$realistic_ranges %||% get_realistic_property_ranges()
  numeric_cols <- intersect(names(ranges), names(simulation_data))

  compute_quality <- function(subset_data) {
    if (nrow(subset_data) == 0 || length(numeric_cols) == 0) return(NA_real_)
    in_range <- vapply(numeric_cols, function(col) {
      r <- ranges[[col]]
      mean(subset_data[[col]] >= r$min & subset_data[[col]] <= r$max, na.rm = TRUE)
    }, numeric(1))
    mean(in_range, na.rm = TRUE)
  }

  if (!("hzdept_r" %in% names(simulation_data))) {
    return(list(surface_horizon_quality = NA_real_, subsurface_quality = NA_real_,
                transition_quality = NA_real_))
  }

  surface <- simulation_data[!is.na(simulation_data$hzdept_r) & simulation_data$hzdept_r <= 30, ]
  subsurface <- simulation_data[!is.na(simulation_data$hzdept_r) & simulation_data$hzdept_r > 30, ]

  transition_quality <- NA_real_
  if (length(numeric_cols) > 0) {
    prop <- numeric_cols[1]
    depth_means <- tapply(simulation_data[[prop]], simulation_data$hzdept_r, mean, na.rm = TRUE)
    depth_means <- depth_means[order(as.numeric(names(depth_means)))]
    depth_means <- depth_means[is.finite(depth_means)]

    if (length(depth_means) >= 2) {
      value_range <- diff(range(depth_means))
      max_jump <- max(abs(diff(depth_means)))
      transition_quality <- if (value_range > 0) max(0, 1 - max_jump / value_range) else 1
    }
  }

  list(
    surface_horizon_quality = compute_quality(surface),
    subsurface_quality = compute_quality(subsurface),
    transition_quality = transition_quality
  )
}

#' Assess Pedological Relationships
#'
#' Checks a small set of cheap, well-established pedological relationship
#' directions when the relevant properties are present: clay content
#' positively associated with CEC, and organic matter negatively associated
#' with depth.
#'
#' @param simulation_data Simulation data.
#' @param criteria Unused (no per-relationship configuration currently
#'   defined upstream).
#' @return List with `texture_relationship_quality` (fraction of rows with
#'   sand+clay+silt within 5 of 100), `chemical_relationship_quality`
#'   (rescaled clay-CEC correlation), `physical_relationship_quality`
#'   (rescaled OM-depth correlation) - each `NA` when the needed properties
#'   are absent.
assess_pedological_relationships <- function(simulation_data, criteria) {
  safe_cor <- function(x, y) {
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 3 || stats::sd(x[ok]) == 0 || stats::sd(y[ok]) == 0) return(NA_real_)
    stats::cor(x[ok], y[ok])
  }

  texture_quality <- NA_real_
  texture_cols <- intersect(c("sandtotal", "claytotal", "silttotal"), names(simulation_data))
  if (length(texture_cols) == 3) {
    texture_sum <- rowSums(simulation_data[texture_cols], na.rm = TRUE)
    texture_quality <- mean(abs(texture_sum - 100) <= 5, na.rm = TRUE)
  }

  chemical_quality <- NA_real_
  if (all(c("claytotal", "cec7") %in% names(simulation_data))) {
    r <- safe_cor(simulation_data$claytotal, simulation_data$cec7)
    if (is.finite(r)) chemical_quality <- (r + 1) / 2
  }

  physical_quality <- NA_real_
  if (all(c("om", "hzdept_r") %in% names(simulation_data))) {
    r <- safe_cor(simulation_data$om, simulation_data$hzdept_r)
    if (is.finite(r)) physical_quality <- (1 - r) / 2
  }

  list(
    texture_relationship_quality = texture_quality,
    chemical_relationship_quality = chemical_quality,
    physical_relationship_quality = physical_quality
  )
}

#' Validate Simulation Depth Trends
#'
#' Compares depth-binned means of `simulation_data` against `original_data`
#' for shared numeric properties, via the correlation between simulated and
#' original per-depth means.
#'
#' @param simulation_data Simulation data.
#' @param original_data Original data for comparison.
#' @param criteria Unused.
#' @return List with `depth_trend_realism` (mean rescaled correlation across
#'   properties), `trend_violations` (count of properties with correlation
#'   below 0.5), `overall_trend_quality`.
validate_simulation_depth_trends <- function(simulation_data, original_data, criteria) {
  if (is.null(original_data) || !("hzdept_r" %in% names(simulation_data)) ||
      !("hzdept_r" %in% names(original_data))) {
    return(list(depth_trend_realism = NA_real_, trend_violations = NA_integer_,
                overall_trend_quality = "unknown"))
  }

  numeric_cols <- intersect(
    names(simulation_data)[vapply(simulation_data, is.numeric, logical(1))],
    names(original_data)[vapply(original_data, is.numeric, logical(1))]
  )
  numeric_cols <- setdiff(numeric_cols, c("simulation_number", "hzdept_r", "hzdepb_r"))

  trend_scores <- c()
  violations <- 0

  for (col in numeric_cols) {
    sim_means <- tapply(simulation_data[[col]], simulation_data$hzdept_r, mean, na.rm = TRUE)
    orig_means <- tapply(original_data[[col]], original_data$hzdept_r, mean, na.rm = TRUE)

    common_depths <- intersect(names(sim_means), names(orig_means))
    if (length(common_depths) >= 3) {
      sim_vals <- sim_means[common_depths]
      orig_vals <- orig_means[common_depths]
      ok <- is.finite(sim_vals) & is.finite(orig_vals)

      if (sum(ok) >= 3 && stats::sd(sim_vals[ok]) > 0 && stats::sd(orig_vals[ok]) > 0) {
        r <- stats::cor(sim_vals[ok], orig_vals[ok])
        trend_scores <- c(trend_scores, (r + 1) / 2)
        if (r < 0.5) violations <- violations + 1
      }
    }
  }

  realism <- if (length(trend_scores) > 0) mean(trend_scores) else NA_real_

  list(
    depth_trend_realism = realism,
    trend_violations = violations,
    overall_trend_quality = if (is.na(realism)) "unknown" else if (realism >= 0.7) "good" else "poor"
  )
}

# Enhanced calculation functions

#' Calculate Workflow Performance
#'
#' Derives real performance metrics from `components` (fraction of expected
#' components that were successfully extracted, i.e. non-`NULL`) and
#' `validation_results` (fraction of `*_validation` steps that neither were
#' skipped nor failed).
#'
#' @param components Extracted workflow components (named list, possibly
#'   with `NULL` entries for components that could not be extracted).
#' @param validation_results Accumulated validation results, with one
#'   `<step>_validation` entry per pipeline step.
#' @return List with `processing_efficiency`, `validation_coverage`,
#'   `overall_performance`.
calculate_workflow_performance <- function(components, validation_results) {
  n_components <- sum(!vapply(components, is.null, logical(1)))
  expected_components <- length(components)
  processing_efficiency <- if (expected_components > 0) n_components / expected_components else NA_real_

  step_results <- validation_results[grepl("_validation$", names(validation_results))]
  n_steps <- length(step_results)
  n_covered <- sum(vapply(step_results, function(x) {
    !is.null(x) && !isTRUE(x$validation_skipped) && !isTRUE(x$validation_failed)
  }, logical(1)))
  validation_coverage <- if (n_steps > 0) n_covered / n_steps else NA_real_

  overall <- if (is.finite(processing_efficiency) && is.finite(validation_coverage)) {
    mean(c(processing_efficiency, validation_coverage))
  } else {
    NA_real_
  }

  list(
    processing_efficiency = processing_efficiency,
    validation_coverage = validation_coverage,
    overall_performance = if (is.na(overall)) "unknown" else if (overall >= 0.7) "good" else "poor"
  )
}

calculate_monte_carlo_score <- function(mc_validation) {
  # Enhanced Monte Carlo scoring  safe operations
  scores <- c()

  if (!is.null(mc_validation$convergence_assessment$converged) && mc_validation$convergence_assessment$converged) {
    scores <- c(scores, 1.0)
  } else {
    scores <- c(scores, 0.5)
  }

  if (!is.null(mc_validation$coverage_assessment$overall_coverage)) {
    scores <- c(scores, mc_validation$coverage_assessment$overall_coverage / 100)
  }

  if (!is.null(mc_validation$distribution_fidelity$overall_fidelity)) {
    scores <- c(scores, mc_validation$distribution_fidelity$overall_fidelity)
  }

  return(mean(scores, na.rm = TRUE))
}

calculate_correlation_score <- function(corr_validation) {
  scores <- c()

  if (!is.null(corr_validation$preservation_assessment$overall_preservation$correlation_rmse)) {
    preservation_score <- max(0, 1 - corr_validation$preservation_assessment$overall_preservation$correlation_rmse * 10)
    scores <- c(scores, preservation_score)
  }

  if (!is.null(corr_validation$matrix_quality$overall_quality)) {
    scores <- c(scores, corr_validation$matrix_quality$overall_quality)
  }

  return(mean(scores, na.rm = TRUE))
}

calculate_gp_model_score <- function(gp_validation) {
  scores <- c()

  if (!is.null(gp_validation$model_performance$overall_performance$mean_r_squared)) {
    scores <- c(scores, gp_validation$model_performance$overall_performance$mean_r_squared)
  }

  if (!is.null(gp_validation$depth_trend_realism$overall_realism_score)) {
    scores <- c(scores, gp_validation$depth_trend_realism$overall_realism_score)
  }

  return(mean(scores, na.rm = TRUE))
}

calculate_soil_science_score <- function(soil_validation) {
  scores <- c()

  if (!is.null(soil_validation$property_constraints$overall_compliance)) {
    scores <- c(scores, soil_validation$property_constraints$overall_compliance)
  }

  if (!is.null(soil_validation$pedological_relationships$overall_adherence)) {
    scores <- c(scores, soil_validation$pedological_relationships$overall_adherence)
  }

  return(mean(scores, na.rm = TRUE))
}

# Enhanced reporting functions
create_report_content <- function(validation_results, include_plots) {
  list(
    executive_summary = create_executive_summary(validation_results),
    detailed_results = create_detailed_results(validation_results),
    plots = if (include_plots) validation_results$diagnostic_plots else NULL,
    recommendations = validation_results$recommendations,
    metadata = validation_results$validation_metadata
  )
}

#' Render Report Content as Markdown Lines
#'
#' Dependency-free renderer that turns the nested list produced by
#' `create_report_content()` into plain markdown lines. Shared by
#' [generate_markdown_report()], [generate_html_report()], and
#' [generate_pdf_report()]'s dependency-free fallback paths.
#'
#' @param report_content Nested list report content.
#' @return Character vector of markdown lines.
render_report_as_markdown_lines <- function(report_content) {
  lines <- c("# Soil Simulation Validation Report", "")
  for (nm in names(report_content)) {
    lines <- c(lines,
              paste0("## ", gsub("_", " ", nm)),
              "",
              render_value_as_lines(report_content[[nm]]),
              "")
  }
  lines
}

#' Render an Arbitrary Value as Indented Markdown Bullet Lines
#'
#' Recursively renders lists as nested bullets and atomic vectors as a
#' single formatted bullet; non-atomic, non-list values (e.g. a ggplot
#' object) are rendered by class name only, deliberately avoiding
#' `print()`/`format()` on arbitrary objects (which could have side effects).
#'
#' @param x Value to render.
#' @param indent Current indentation depth.
#' @return Character vector of markdown lines.
render_value_as_lines <- function(x, indent = 0) {
  prefix <- strrep("  ", indent)

  if (is.null(x)) {
    return(paste0(prefix, "- (none)"))
  }

  # A "plain" list (no S3 class of its own, e.g. not a data.frame or some
  # other classed object masquerading as a list under the hood) recurses as
  # nested bullets; anything with its own class falls through to the
  # class-name-only branch below instead.
  if (is.list(x) && (is.null(attr(x, "class")) || identical(class(x), "list"))) {
    if (length(x) == 0) return(paste0(prefix, "- (empty)"))
    nm <- names(x)
    if (is.null(nm) || any(nm == "")) nm <- paste0("[", seq_along(x), "]")

    out <- character(0)
    for (i in seq_along(x)) {
      out <- c(out, paste0(prefix, "- ", nm[i], ":"), render_value_as_lines(x[[i]], indent + 1))
    }
    return(out)
  }

  if (is.atomic(x)) {
    return(paste0(prefix, "- ", paste(format(x), collapse = ", ")))
  }

  paste0(prefix, "- <", paste(class(x), collapse = "/"), ">")
}

#' Escape Text for Inclusion in HTML
#'
#' @param x Character vector.
#' @return Character vector with `&`, `<`, `>` escaped.
escape_html_text <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x
}

#' Generate HTML Validation Report
#'
#' Writes `report_content` to `output_file` as HTML. Uses
#' `rmarkdown::render()` for richer output when `rmarkdown` is installed and
#' `pandoc` is available (both optional, `Suggests`-only dependencies);
#' otherwise falls back to a dependency-free `<pre>`-wrapped rendering that
#' always works.
#'
#' @param report_content Report content from `create_report_content()`.
#' @param output_file Output file path (`.html`).
#' @return `TRUE` on success, `FALSE` on failure (logged as a warning).
generate_html_report <- function(report_content, output_file) {
  tryCatch({
    md_lines <- render_report_as_markdown_lines(report_content)

    if (requireNamespace("rmarkdown", quietly = TRUE) && nzchar(Sys.which("pandoc"))) {
      tmp_md <- tempfile(fileext = ".md")
      writeLines(md_lines, tmp_md)
      rmarkdown::render(tmp_md, output_format = "html_document",
                        output_file = basename(output_file),
                        output_dir = dirname(output_file), quiet = TRUE)
    } else {
      html_lines <- c("<!doctype html><html><body><pre>",
                      escape_html_text(md_lines),
                      "</pre></body></html>")
      writeLines(html_lines, output_file)
    }
    TRUE
  }, error = function(e) {
    handle_workflow_error(e, "HTML report generation", "warn")
    FALSE
  })
}

#' Generate Markdown Validation Report
#'
#' Writes `report_content` to `output_file` as a plain markdown file
#' (dependency-free).
#'
#' @param report_content Report content from `create_report_content()`.
#' @param output_file Output file path (`.md`).
#' @return `TRUE` on success, `FALSE` on failure (logged as a warning).
generate_markdown_report <- function(report_content, output_file) {
  tryCatch({
    writeLines(render_report_as_markdown_lines(report_content), output_file)
    TRUE
  }, error = function(e) {
    handle_workflow_error(e, "Markdown report generation", "warn")
    FALSE
  })
}

#' Generate PDF Validation Report
#'
#' Writes `report_content` to `output_file` as a PDF. Uses
#' `rmarkdown::render()` for richer output when `rmarkdown`/`tinytex` are
#' installed, `pandoc` is available, and a working TeX distribution is
#' detected (all optional, `Suggests`-only dependencies); otherwise falls
#' back to a dependency-free plain-text rendering via [grDevices::pdf()]
#' that always works.
#'
#' @param report_content Report content from `create_report_content()`.
#' @param output_file Output file path (`.pdf`).
#' @return `TRUE` on success, `FALSE` on failure (logged as a warning).
generate_pdf_report <- function(report_content, output_file) {
  tryCatch({
    md_lines <- render_report_as_markdown_lines(report_content)

    can_use_rmarkdown <- requireNamespace("rmarkdown", quietly = TRUE) &&
      requireNamespace("tinytex", quietly = TRUE) &&
      nzchar(Sys.which("pandoc")) &&
      isTRUE(tryCatch(tinytex::is_tinytex(), error = function(e) FALSE))

    if (can_use_rmarkdown) {
      tmp_md <- tempfile(fileext = ".md")
      writeLines(md_lines, tmp_md)
      rmarkdown::render(tmp_md, output_format = "pdf_document",
                        output_file = basename(output_file),
                        output_dir = dirname(output_file), quiet = TRUE)
    } else {
      grDevices::pdf(output_file)
      on.exit(grDevices::dev.off(), add = TRUE)

      lines_per_page <- 45
      n_pages <- max(1, ceiling(length(md_lines) / lines_per_page))

      for (p in seq_len(n_pages)) {
        graphics::plot.new()
        page_lines <- md_lines[(((p - 1) * lines_per_page) + 1):min(p * lines_per_page, length(md_lines))]
        graphics::text(0, 1, paste(page_lines, collapse = "\n"),
                       adj = c(0, 1), family = "mono", cex = 0.7)
      }
    }
    TRUE
  }, error = function(e) {
    handle_workflow_error(e, "PDF report generation", "warn")
    FALSE
  })
}

#' Shared minimal ggplot theme for diagnostic plots
#'
#' No ggplot styling convention existed anywhere in modules/ prior to this
#' (grepped: zero ggplot()/geom_*() calls anywhere despite ggplot2 being
#' loaded in this file and mod06_gp_modeling.R) - this establishes one
#' rather than improvising ad hoc per plot below.
theme_soil_diagnostics <- function() {
  ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom"
    )
}

#' Save a ggplot object to output_dir if one was provided
#' @param plot_obj A ggplot object to save.
#' @param filename File name (not path) to save the plot under.
#' @param output_dir Directory to save into, or `NULL` to skip saving.
#' @return The ggplot object itself (output_dir NULL or save failed) or the
#'   saved file path (output_dir provided and save succeeded).
save_diagnostic_plot <- function(plot_obj, filename, output_dir) {
  if (is.null(output_dir)) {
    return(plot_obj)
  }
  tryCatch({
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    file_path <- file.path(output_dir, filename)
    ggplot2::ggsave(file_path, plot = plot_obj, width = 8, height = 6)
    file_path
  }, error = function(e) {
    log_message("WARN", paste("Failed to save plot", filename, ":", conditionMessage(e)), category = "Validation")
    plot_obj
  })
}

#' Generate diagnostic plots across all validated workflow components
#'
#' Builds real plots from `components` (mod01-07 outputs) and
#' `validation_results` (this file's own already-computed assessments)
#' rather than returning empty placeholders. Each of the 5 categories is
#' independently tryCatch-wrapped so a missing/malformed component degrades
#' that one category to an empty list rather than failing the whole call -
#' matching this file's existing per-step error-tolerance convention.
#'
#' @param components From extract_and_validate_components(): may contain
#'   simulation_data, monte_carlo_results, final_data, gp_models,
#'   correlation_matrices, training_data.
#' @param validation_results The accumulating validation_results list built
#'   by execute_validation_pipeline() (soil_science_validation,
#'   overall_assessment, etc.).
#' @param output_dir If non-NULL, each plot is additionally saved there
#'   (nothing downstream currently consumes these plot objects otherwise -
#'   generate_html_report()/generate_markdown_report()/generate_pdf_report()
#'   are themselves still-unimplemented stubs).
#' @return List with monte_carlo_plots, correlation_plots, gp_model_plots,
#'   soil_science_plots, summary_plots (each a named list of plot objects or
#'   file paths, empty if that category's data wasn't available), and
#'   plot_generation_status.
generate_comprehensive_diagnostics <- function(components, validation_results, output_dir) {

  # ggplot2 is Suggests-only (not a hard package dependency): every plot
  # category below that renders via ggplot2 checks this flag first and
  # degrades to an empty list when it's unavailable, rather than only
  # relying on the surrounding tryCatch to catch the resulting
  # "there is no package called 'ggplot2'" error.
  ggplot2_available <- requireNamespace("ggplot2", quietly = TRUE)

  get_simulation_data <- function() {
    if (is.data.frame(components$simulation_data)) {
      return(components$simulation_data)
    }
    if (is.data.frame(components$monte_carlo_results$simulation_data)) {
      return(components$monte_carlo_results$simulation_data)
    }
    NULL
  }

  numeric_property_columns <- function(sim_data) {
    exclude_cols <- c("cokey", "simulation_number", "hzdept_r", "hzdepb_r")
    setdiff(names(sim_data)[vapply(sim_data, is.numeric, logical(1))], exclude_cols)
  }

  # 1. Monte Carlo simulated-property distribution plots
  monte_carlo_plots <- tryCatch({
    sim_data <- get_simulation_data()
    props <- if (!is.null(sim_data)) numeric_property_columns(sim_data) else character(0)

    if (!ggplot2_available || is.null(sim_data) || length(props) == 0) {
      log_message("WARN", "No simulation_data (or ggplot2 unavailable) for Monte Carlo diagnostic plots", category = "Validation")
      list()
    } else {
      long_data <- tidyr::pivot_longer(sim_data[props], cols = dplyr::everything(),
                                        names_to = "property", values_to = "value")
      p <- ggplot2::ggplot(long_data, ggplot2::aes(x = value)) +
        ggplot2::geom_histogram(bins = 30, fill = "#4c72b0", color = "white", na.rm = TRUE) +
        ggplot2::facet_wrap(~property, scales = "free") +
        ggplot2::labs(title = "Monte Carlo Simulated Property Distributions", x = NULL, y = "Count") +
        theme_soil_diagnostics()
      list(property_distributions = save_diagnostic_plot(p, "monte_carlo_distributions.png", output_dir))
    }
  }, error = function(e) {
    log_message("WARN", paste("Monte Carlo diagnostic plot generation failed:", conditionMessage(e)), category = "Validation")
    list()
  })

  # 2. Correlation structure before/after (corrplot draws to a base-graphics
  # device rather than returning an object, so this always renders to a PNG
  # file - output_dir if given, otherwise a temp file - and returns the path.
  correlation_plots <- tryCatch({
    sim_data <- get_simulation_data()
    orig_cor <- components$correlation_matrices$global_correlation_matrix

    if (is.null(sim_data) || is.null(orig_cor) || !requireNamespace("corrplot", quietly = TRUE)) {
      log_message("WARN", "Missing simulation_data/correlation_matrices for correlation diagnostic plot", category = "Validation")
      list()
    } else {
      props <- numeric_property_columns(sim_data)
      sim_cor <- safe_correlation(sim_data[props], method = "pearson", handle_constant = "warn")
      common_props <- intersect(rownames(orig_cor), rownames(sim_cor))

      if (length(common_props) < 2) {
        list()
      } else {
        orig_subset <- orig_cor[common_props, common_props]
        sim_subset <- sim_cor[common_props, common_props]

        file_path <- file.path(if (!is.null(output_dir)) output_dir else tempdir(), "correlation_before_after.png")
        if (!is.null(output_dir) && !dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

        grDevices::png(file_path, width = 1000, height = 500)
        graphics::par(mfrow = c(1, 2))
        corrplot::corrplot(orig_subset, method = "color", title = "Original", mar = c(0, 0, 2, 0))
        corrplot::corrplot(sim_subset, method = "color", title = "Simulated", mar = c(0, 0, 2, 0))
        grDevices::dev.off()

        list(before_after = file_path)
      }
    }
  }, error = function(e) {
    log_message("WARN", paste("Correlation diagnostic plot generation failed:", conditionMessage(e)), category = "Validation")
    list()
  })

  # 3. GP depth-trend plots (recomputed from components$gp_models the same
  # way assess_depth_trend_realism() does, since that function's own
  # trend_predictions doesn't retain the raw predicted values - see
  # assess_trend_realism(), which discards `predictions` entirely).
  gp_model_plots <- tryCatch({
    gp_models <- components$gp_models
    if (!ggplot2_available || is.null(gp_models)) {
      log_message("WARN", "No gp_models (or ggplot2 unavailable) for GP diagnostic plots", category = "Validation")
      list()
    } else {
      test_depths <- seq(0, 200, by = 5)
      properties <- setdiff(names(gp_models), "model_summary")

      rows <- list()
      for (prop in properties) {
        prop_models <- gp_models[[prop]]
        if (is.null(prop_models$models)) next
        for (group in names(prop_models$models)) {
          group_model <- prop_models$models[[group]]
          if (is.null(group_model)) next
          predicted <- tryCatch(
            predict_gp_depth_trends(group_model, test_depths),
            error = function(e) rep(NA_real_, length(test_depths))
          )
          rows[[paste(prop, group)]] <- data.frame(
            property = prop, group = group, depth = test_depths, predicted = predicted
          )
        }
      }

      if (length(rows) == 0) {
        list()
      } else {
        trend_data <- do.call(rbind, rows)
        p <- ggplot2::ggplot(trend_data, ggplot2::aes(x = depth, y = predicted, color = group)) +
          ggplot2::geom_line(na.rm = TRUE) +
          ggplot2::facet_wrap(~property, scales = "free_y") +
          ggplot2::labs(title = "GP Depth-Trend Predictions by Property/Group", x = "Depth (cm)", y = "Predicted value") +
          theme_soil_diagnostics()
        list(depth_trends = save_diagnostic_plot(p, "gp_depth_trends.png", output_dir))
      }
    }
  }, error = function(e) {
    log_message("WARN", paste("GP diagnostic plot generation failed:", conditionMessage(e)), category = "Validation")
    list()
  })

  # 4. Soil-science range-violation plot. range_violations is the one part
  # of assess_property_constraints()'s output that's genuinely data-derived
  # (cross_property_violations/distribution_anomalies are themselves still
  # hardcoded stubs elsewhere in this file - not plotted here since they
  # aren't real data).
  soil_science_plots <- tryCatch({
    range_violations <- validation_results$soil_science_validation$property_constraints$range_violations
    if (!ggplot2_available || is.null(range_violations) || length(range_violations) == 0) {
      log_message("WARN", "No range_violations (or ggplot2 unavailable) for soil science diagnostic plots", category = "Validation")
      list()
    } else {
      violation_data <- data.frame(
        property = names(range_violations),
        pct_violations = vapply(range_violations, function(v) v$pct_violations %||% NA_real_, numeric(1))
      )
      p <- ggplot2::ggplot(violation_data, ggplot2::aes(x = stats::reorder(property, pct_violations), y = pct_violations)) +
        ggplot2::geom_col(fill = "#c44e52") +
        ggplot2::coord_flip() +
        ggplot2::labs(title = "Property Range Violations", x = NULL, y = "% of values out of range") +
        theme_soil_diagnostics()
      list(range_violations = save_diagnostic_plot(p, "soil_science_range_violations.png", output_dir))
    }
  }, error = function(e) {
    log_message("WARN", paste("Soil science diagnostic plot generation failed:", conditionMessage(e)), category = "Validation")
    list()
  })

  # 5. Overall quality-score summary plot
  summary_plots <- tryCatch({
    component_scores <- validation_results$overall_assessment$component_scores
    if (!ggplot2_available || is.null(component_scores) || length(component_scores) == 0) {
      log_message("WARN", "No component_scores (or ggplot2 unavailable) for summary diagnostic plot", category = "Validation")
      list()
    } else {
      score_data <- data.frame(
        component = names(component_scores),
        score = vapply(component_scores, function(s) s %||% NA_real_, numeric(1))
      )
      overall_score <- validation_results$overall_assessment$quality_score
      overall_grade <- validation_results$overall_assessment$quality_grade
      subtitle <- if (!is.null(overall_score)) {
        sprintf("Overall: %.2f (%s)", overall_score, overall_grade %||% "")
      } else {
        NULL
      }
      p <- ggplot2::ggplot(score_data, ggplot2::aes(x = component, y = score)) +
        ggplot2::geom_col(fill = "#55a868") +
        ggplot2::ylim(0, 1) +
        ggplot2::labs(title = "Workflow Quality Scores by Component", subtitle = subtitle, x = NULL, y = "Score (0-1)") +
        theme_soil_diagnostics()
      list(quality_scores = save_diagnostic_plot(p, "quality_scores_summary.png", output_dir))
    }
  }, error = function(e) {
    log_message("WARN", paste("Summary diagnostic plot generation failed:", conditionMessage(e)), category = "Validation")
    list()
  })

  list(
    monte_carlo_plots = monte_carlo_plots,
    correlation_plots = correlation_plots,
    gp_model_plots = gp_model_plots,
    soil_science_plots = soil_science_plots,
    summary_plots = summary_plots,
    plot_generation_status = "complete"
  )
}

generate_workflow_recommendations <- function(validation_results) {
  recommendations <- character(0)

  # Generate recommendations based on validation results
  if (!is.null(validation_results$overall_assessment$quality_score) &&
      validation_results$overall_assessment$quality_score < 0.8) {
    recommendations <- c(recommendations, "Consider improving overall workflow quality")
  }

  if (length(validation_results$overall_assessment$critical_issues) > 0) {
    recommendations <- c(recommendations, "Address critical issues identified in validation")
  }

  return(recommendations)
}

# Enhanced utility functions
get_default_quality_weights <- function() {
  list(
    monte_carlo = 0.25,
    correlation = 0.25,
    gp_models = 0.25,
    soil_science = 0.25
  )
}

get_realistic_property_ranges <- function() {
  # NOTE: the original source called get_predefined_properties("laboratory")
  # here via `|> {list(...)}`, but the block never referenced `.` (the
  # piped-in value) - a no-op pipe. The call was dropped; this is
  # behavior-preserving since get_predefined_properties() is a pure lookup
  # with no side effects and its result was never used.
  list(
    clay_total = list(min = 0, max = 100),
    sand_total = list(min = 0, max = 100),
    silt_total = list(min = 0, max = 100),
    ph = list(min = 2.5, max = 11.0),
    bulk_density = list(min = 0.3, max = 3.0),
    water_retention_third_bar = list(min = 0, max = 70),
    water_retention_15_bar = list(min = 0, max = 50),
    om = list(min = 0, max = 50),
    cec = list(min = 0, max = 200),
    rfv = list(min = 0, max = 95)
  )
}

get_realistic_gradients <- function() {
  list(
    clay_total = list(max_gradient = 5),
    sand_total = list(max_gradient = 5),
    ph = list(max_gradient = 0.5),
    bulk_density = list(max_gradient = 0.2)
  )
}

create_executive_summary <- function(validation_results) {
  list(
    overall_score = validation_results$overall_assessment$quality_score,
    workflow_status = validation_results$overall_assessment$workflow_status,
    key_findings = "Enhanced validation completed",
    critical_issues = validation_results$overall_assessment$critical_issues
  )
}

create_detailed_results <- function(validation_results) {
  list(
    component_scores = validation_results$overall_assessment$component_scores,
    validation_details = "Detailed results available",
    methodology = "Enhanced validation  utilities"
  )
}

#' Perform GP Cross-Validation
#'
#' Runs real k-fold cross-validation (reusing the shared [k_fold_gp_cv()]
#' helper from `gp-modeling.R`) for each property in `gp_models` against the
#' matching column of `training_data`, aggregating mean CV RMSE and an
#' R-squared derived from it.
#'
#' @param gp_models GP models from `build_stratified_gp_models()`.
#' @param training_data Original training data with `hzdept_r` and one
#'   column per property in `gp_models`.
#' @param criteria List, optionally with `n_folds` (default 5).
#' @return List with `cv_rmse`, `cv_r_squared`, `cv_quality`, `cv_method`.
perform_gp_cross_validation <- function(gp_models, training_data, criteria) {
  n_folds <- criteria$n_folds %||% 5
  properties <- names(gp_models)[names(gp_models) != "model_summary"]

  fold_rmses <- c()
  fold_r2 <- c()

  for (prop in properties) {
    if (is.null(training_data) || !(prop %in% names(training_data)) ||
        !("hzdept_r" %in% names(training_data))) next

    depths <- training_data$hzdept_r
    values <- training_data[[prop]]
    ok <- is.finite(depths) & is.finite(values)

    min_needed <- 2 * max(2, min(n_folds, sum(ok)))
    if (sum(ok) < min_needed) next

    depth_range <- range(depths[ok])
    if (diff(depth_range) <= 0) next
    scaled_depths <- (depths[ok] - depth_range[1]) / diff(depth_range)

    cv <- k_fold_gp_cv(scaled_depths, values[ok], n_folds,
                       list(default = list(type = "exponential", power = 1.95)))
    rmse <- cv$mean_rmse_by_candidate$default

    if (is.finite(rmse)) {
      fold_rmses <- c(fold_rmses, rmse)
      total_var <- stats::var(values[ok])
      if (is.finite(total_var) && total_var > 0) {
        fold_r2 <- c(fold_r2, max(0, 1 - (rmse^2 / total_var)))
      }
    }
  }

  cv_rmse <- if (length(fold_rmses) > 0) mean(fold_rmses) else NA_real_
  cv_r_squared <- if (length(fold_r2) > 0) mean(fold_r2) else NA_real_

  list(
    cv_rmse = cv_rmse,
    cv_r_squared = cv_r_squared,
    cv_quality = if (is.na(cv_r_squared)) "unknown" else if (cv_r_squared >= 0.6) "good" else "poor",
    cv_method = "k_fold_gp_cv"
  )
}

