# Generated from utils.R (P1b split): file I/O & data management


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
#' @section Usage note:
#' Fully implemented and exported, but not currently called from anywhere
#' else in the package outside a \code{@seealso} doc link (confirmed by
#' source grep) - standalone public API for callers who need it, not
#' dead/broken code. See also \code{\link{write_soil_data}}.
#' @keywords internal
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
#' @section Usage note:
#' Fully implemented and exported, but not currently called from anywhere
#' else in the package outside a \code{@seealso} doc link (confirmed by
#' source grep) - standalone public API for callers who need it, not
#' dead/broken code. See also \code{\link{read_soil_data}}.
#' @keywords internal
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
#' @section Usage note:
#' Fully implemented and exported, but not currently called from anywhere
#' else in the package (confirmed by source grep) - standalone public API
#' for callers who need standalone/versioned backups (as opposed to
#' \code{\link{write_soil_data}}'s single-backup-on-overwrite behavior), not
#' dead/broken code.
#' @keywords internal
backup_data <- function(source_path, backup_dir = NULL, max_backups = 5) {

  if (!file.exists(source_path)) {
    stop("Source file not found: ", source_path)
  }

  if (is.null(backup_dir)) {
    backup_dir <- dirname(source_path)
  }

  dir.create(backup_dir, recursive = TRUE, showWarnings = FALSE)

  # Create timestamped backup filename (shared with write_soil_data())
  base_name <- tools::file_path_sans_ext(basename(source_path))
  backup_path <- create_backup_filename(source_path, backup_dir)

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
#' @keywords internal
load_configuration <- function(config_path, config_type = "auto", validate_config = TRUE,
                               verbose = getOption("ssurgo.verbose", FALSE)) {
  .old_log_cfg <- set_verbose_logging(verbose)
  on.exit(options(soil_workflow_log_config = .old_log_cfg), add = TRUE)

  if (!file.exists(config_path)) {
    log_message("WARN", paste("Configuration file not found:", config_path, ". Using defaults."), category = "Config")
    return(default_config())
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
  default_config <- default_config()
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

