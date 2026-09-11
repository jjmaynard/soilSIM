# Guardrail for the naming conventions in CONTRIBUTING.md.
#
# `approved_prefix_regex` encodes the accepted verb/noun lexicon. `legacy_exceptions`
# is the shrinking list of pre-0.2.0 names that don't fit it yet and haven't been
# renamed (each needs a deprecation shim when it is). The test fails if a NEW
# offender appears or if a name leaves the exception list without being removed
# from the package.

approved_prefix_regex <- paste0(
  "^(fetch|read|write|fit|quantile|simulate|fuse|estimate|compute|analyze|",
  "validate|diagnose|default|convert|prepare|adjust|apply|infill|process|",
  "build|lookup|rasterize|resolve|classify|remarginalize|extract|update|",
  "attach|align|correct|impute|integrate|learn|match|merge|normalize|parse|",
  "predict|restore|run|select|setup|slice|standardize|summarize|track|",
  "new|cache|load|save|detect|clean|download|remove|sample|query|ensure|",
  "identify|create|generate|compare|evaluate|handle|backup|add|assess|",
  "aggregate|is|has)_",
  "|_raster$",
  "|^(available_properties|predefined_properties|distributions_for_properties|",
  "rfv_range_category|group_members|ilr_forward|ilr_inverse|van_genuchten|",
  "tri_dist|property_contextual_ranges|CACHE_TTL_SECONDS|SSURGO_SIM_PROPERTY_COLUMNS)$"
)

# names defined only as deprecation shims are exempt (they intentionally keep the old name)
deprecated_shim_names <- function() {
  f <- testthat::test_path("..", "..", "R", "deprecated.R")
  if (!file.exists(f)) return(character(0))
  txt <- readLines(f, warn = FALSE)
  unique(regmatches(txt, regexpr("^[a-zA-Z0-9_.]+(?= <- function)", txt, perl = TRUE)))
}

# Residual legacy names (pre-0.2.0) that still need a rename + shim. Only shrink this.
# (D1's export-tiering pass demoted most of the original 23 to @keywords internal,
# so they no longer appear in getNamespaceExports() and dropped out of this list on
# their own; these 8 are still public and still need the D2 rename.)
legacy_exceptions <- c(
  "beta_to_moments", "gamma_to_moments", "moments_to_beta", "moments_to_gamma",
  "lognormal_to_normal_params", "normal_to_lognormal_params",
  "unwrap_nested_rasters", "remarginalized_awc"
)

offending_exports <- function() {
  ns <- asNamespace("soilSIM")
  exports <- getNamespaceExports("soilSIM")
  fns <- exports[vapply(
    exports,
    function(x) is.function(get0(x, envir = ns, inherits = FALSE)),
    logical(1)
  )]
  setdiff(fns[!grepl(approved_prefix_regex, fns)], deprecated_shim_names())
}

test_that("no new naming-convention offenders", {
  offenders <- offending_exports()
  new_offenders <- setdiff(offenders, legacy_exceptions)
  expect_equal(new_offenders, character(0))
})

test_that("legacy exception list only shrinks", {
  offenders <- offending_exports()
  stale <- setdiff(legacy_exceptions, offenders)
  expect_equal(
    stale, character(0),
    info = "these names are no longer offenders; remove them from legacy_exceptions"
  )
})

test_that("public (non-internal) export count trend", {
  skip("hard cap enabled once export demotion (P5) lands")
  expect_lte(length(getNamespaceExports("soilSIM")), 120)
})
