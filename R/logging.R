# Generated from utils.R (P1b split): logging & error handling


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

