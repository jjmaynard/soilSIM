# Make lifecycle deprecation conditions loud during the test run so the shims in
# R/deprecated.R are actually exercised (deprecate_soft() is otherwise silent
# outside the global environment).
withr::local_options(
  lifecycle_verbosity = "warning",
  .local_envir = testthat::teardown_env()
)
