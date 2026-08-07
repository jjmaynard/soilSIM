#' @title Shared `future`/`future.apply` Parallel Dispatch Helper
#' @description A single internal helper standardizing this package's parallel processing on
#'   `future`/`future.apply` (`future::multisession`, which behaves identically on Windows and
#'   Unix), replacing the hand-rolled base-`parallel::` boilerplate (OS branching between
#'   `parallel::makeCluster()`/`parLapply()` on Windows and `parallel::mclapply()` elsewhere,
#'   explicit `library(soilSIM)` on Windows workers, explicit `soil_workflow_log_config`
#'   propagation) that used to be duplicated across `process_cokeys_parallel()`
#'   (`multivariate-adjustment.R`), `run_parallel_simulation()` (`monte-carlo.R`), and
#'   `maybe_adjust_soil_data_depth_trend()` (`ssurgo-simulation.R`).
#' @name parallel_utils
NULL

#' Run a Function Over a List/Vector, in Parallel via `future.apply`, with Sequential Fallback
#'
#' @param X A list or vector to iterate over.
#' @param FUN A function taking one element of `X` as its first argument, plus `...`.
#' @param ... Additional fixed arguments forwarded to `FUN` on every call.
#' @param n_cores Number of parallel workers. `NULL` (default) resolves to
#'   `max(1, parallel::detectCores() - 1)`, then is clamped to `min(n_cores, length(X))` (no point
#'   starting more workers than there is work). A resolved value `<= 1` skips `future::plan()`
#'   setup entirely and just runs `lapply(X, FUN, ...)`.
#' @param future_seed Passed through to `future.apply::future_lapply()`'s `future.seed` argument.
#'   Use `TRUE` whenever `FUN` (or anything it calls) uses R's RNG, to get `future`'s
#'   parallel-safe RNG streams - this is more callers than it first appears: besides the obvious
#'   cases (Monte Carlo draws, profile perturbation), `GPfit::GP_fit()`'s internal hyperparameter
#'   search is a genetic algorithm and genuinely uses RNG too (confirmed empirically via
#'   `future_lapply()`'s "UNRELIABLE VALUE" warning when this was left `FALSE` for a GP-fitting
#'   caller in this package), even though its fitted result is highly stable across runs. Default
#'   `FALSE` matches `future_lapply()`'s own default; only leave it `FALSE` for `FUN` you've
#'   confirmed doesn't touch the RNG anywhere in its call chain.
#' @param op_name A short label used in the log message and as `handle_workflow_error()`'s
#'   context string if dispatch fails.
#' @param sequential_fallback A zero-argument closure invoked if dispatch fails and `catch_errors
#'   = TRUE` (the caller pre-binds its own arguments, e.g. `function() run_sequential_simulation(
#'   simulation_params, ...)`). `NULL` (default) falls back to plain `lapply(X, FUN, ...)` instead.
#'   Ignored when `catch_errors = FALSE`.
#' @param catch_errors If `TRUE` (default), wraps dispatch in `tryCatch()` and applies the
#'   `sequential_fallback` contract above on error. If `FALSE`, dispatch errors propagate to the
#'   caller's own error handling instead (plan setup/restoration via `on.exit()` still happens
#'   regardless) - needed for callers whose own contract is "no sequential fallback, propagate/
#'   convert the error myself" (e.g. `simulate_profile_depths_by_collection_parallel()`, which
#'   returns `NULL` on error via its own outer `tryCatch()`).
#' @return A list of `FUN(x, ...)` results, one per element of `X`, in `X`'s order - the same
#'   shape `parallel::parLapply()`/`mclapply()` already produced, so callers that previously
#'   consumed those results (e.g. via `dplyr::bind_rows()` or manual concatenation) need no change.
#' @keywords internal
run_parallel_lapply <- function(X, FUN, ..., n_cores = NULL, future_seed = FALSE,
                                 op_name = "Parallel operation",
                                 sequential_fallback = NULL, catch_errors = TRUE) {
  if (is.null(n_cores)) {
    n_cores <- max(1, parallel::detectCores() - 1)
  }
  n_cores <- min(n_cores, length(X))

  if (n_cores <= 1) {
    return(lapply(X, FUN, ...))
  }

  log_message("INFO", paste("Running", op_name, "in parallel with", n_cores, "workers"), category = "Parallel")

  # Snapshot the CALLER's plan (not a hardcoded assumption of sequential()) and restore it
  # unconditionally via on.exit() - fixes a real bug in an earlier future-based function in this
  # package, which only reverted to sequential() on its own success path, clobbering a caller's
  # pre-existing plan (or leaving multisession active indefinitely) on error.
  old_plan <- future::plan()
  future::plan(future::multisession, workers = n_cores)
  on.exit(future::plan(old_plan), add = TRUE)

  # future::multisession workers are fresh R processes and do not inherit this process's
  # options() any more than parallel::makeCluster() PSOCK workers did - propagate the (possibly
  # verbose-raised) log config explicitly so log_message() calls inside FUN respect the caller's
  # verbose setting. future.apply has no clusterCall()-equivalent "run once at pool creation" hook,
  # so this is applied redundantly per dispatched element rather than once per worker - cheap and
  # correct either way.
  #
  # No clusterExport()/library(soilSIM)-equivalent step is needed here: future's automatic globals
  # detection (via the globals package) exports FUN's referenced variables and attaches the
  # soilSIM namespace automatically when FUN calls a soilSIM-namespaced function - this only works
  # when soilSIM is actually installed and attached (library(soilSIM)) in the calling session, not
  # merely devtools::load_all()'d.
  current_log_cfg <- getOption("soil_workflow_log_config")
  wrapped_FUN <- function(x, ...) {
    options(soil_workflow_log_config = current_log_cfg)
    FUN(x, ...)
  }

  dispatch <- function() {
    future.apply::future_lapply(X, wrapped_FUN, ..., future.seed = future_seed)
  }

  if (!catch_errors) {
    return(dispatch())
  }

  tryCatch(dispatch(), error = function(e) {
    handle_workflow_error(e, op_name, "warn")
    if (!is.null(sequential_fallback)) {
      sequential_fallback()
    } else {
      lapply(X, FUN, ...)
    }
  })
}
