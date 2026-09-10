# Every renamed export keeps a forwarding shim in R/deprecated.R. This checks a
# representative sample: each shim (a) emits a lifecycle deprecation condition and
# (b) returns exactly what the new name returns for the same inputs.
#
# setup.R sets lifecycle_verbosity = "warning" so deprecate_soft() actually fires.

expect_shim <- function(old_call, new_call) {
  old_call <- rlang::enquo(old_call)
  new_call <- rlang::enquo(new_call)

  # (a) the shim is wired to lifecycle: forcing "error" verbosity makes the
  #     deprecation signal throw a classed condition every call.
  withr::with_options(
    list(lifecycle_verbosity = "error"),
    expect_error(rlang::eval_tidy(old_call), class = "lifecycle_error_deprecated")
  )

  # (b) with deprecation quiet, the shim forwards its result unchanged.
  withr::with_options(list(lifecycle_verbosity = "quiet"), {
    withr::local_seed(1)
    old_val <- rlang::eval_tidy(old_call)
    withr::local_seed(1)
    new_val <- rlang::eval_tidy(new_call)
    expect_equal(old_val, new_val)
  })
}

test_that("config-accessor shims forward and warn", {
  expect_shim(get_default_configuration(), default_config())
  expect_shim(get_monte_carlo_defaults(), default_monte_carlo_config())
  expect_shim(get_statistical_analysis_defaults(), default_statistics_config())
  expect_shim(get_validation_defaults("geographic", "epsg:4326"),
              default_diagnostics_config("geographic", "epsg:4326"))
  expect_shim(get_default_property_config("claytotal"), default_property_config("claytotal"))
  expect_shim(get_predefined_properties("ssurgo"), predefined_properties("ssurgo"))
  expect_shim(get_rfv_range_category(15), rfv_range_category(15))
})

test_that("numeric-helper shims forward and warn", {
  expect_shim(calculate_mode(c(1, 2, 2, 3)), compute_mode(c(1, 2, 2, 3)))
  expect_shim(calculate_complexity_score(50, 1, "POLYGON"),
              compute_complexity_score(50, 1, "POLYGON"))
})

test_that("Bayesian-fusion shims forward and warn", {
  expect_shim(bayes_update_normal_normal(10, 2, 12, 3), fuse_normal_normal(10, 2, 12, 3))
  expect_shim(
    bayes_fuse(list(mu = 10, sigma = 2), list(mu = 12, sigma = 3), family = "normal"),
    fuse_distribution(list(mu = 10, sigma = 2), list(mu = 12, sigma = 3), family = "normal")
  )
})

test_that("Tier-1 rename shims warn (smoke, offline)", {
  # These need heavy fixtures to run end to end; here we only assert the shim
  # wrapper itself is a deprecate_soft forwarder, by inspecting its body.
  for (nm in c("generate_monte_carlo_realizations", "build_stratified_gp_models",
               "run_stage1_fusion", "calculate_aws_df", "validate_complete_workflow",
               "download_and_prepare_ssurgo", "sim_component_comp")) {
    b <- paste(deparse(body(get(nm, envir = asNamespace("soilSIM")))), collapse = " ")
    expect_match(b, "deprecate_soft", info = nm)
  }
})

test_that("every deprecated.R shim is exported and resolves", {
  shim_src <- readLines(testthat::test_path("..", "..", "R", "deprecated.R"))
  shims <- unique(regmatches(
    shim_src, regexpr("^[a-zA-Z0-9_.]+(?= <- function)", shim_src, perl = TRUE)
  ))
  expect_gt(length(shims), 40)
  exports <- getNamespaceExports("soilSIM")
  expect_setequal(intersect(shims, exports), shims)
})
