# Guardrail for the naming conventions documented in CONTRIBUTING.md.
#
# Enforcement is deferred until phase P6 of the naming re-architecture
# (see planning-docs / the approved rename plan). Until then this test is
# skipped but the logic is kept live so `offending_exports()` can be run
# interactively to watch the offender list shrink phase by phase.

approved_prefix_regex <- paste0(
  "^(fetch|read|write|fit|quantile|simulate|fuse|estimate|compute|analyze|",
  "validate|diagnose|default|convert|prepare|adjust|apply|infill|process|",
  "build|lookup|rasterize|resolve|classify|remarginalize|extract|update|",
  "attach|align|correct|impute|integrate|learn|match|merge|normalize|parse|",
  "predict|restore|run|select|setup|slice|standardize|summarize|track|",
  "van_genuchten|has|is)_",
  "|_raster$",
  "|^(available|predefined|distributions_for|rfv_range)_",
  "|^(CACHE_TTL_SECONDS|SSURGO_SIM_PROPERTY_COLUMNS)$"
)

deprecated_alias_names <- function() {
  f <- testthat::test_path("..", "..", "R", "deprecated.R")
  if (!file.exists(f)) return(character(0))
  txt <- readLines(f, warn = FALSE)
  m <- regmatches(txt, regexpr("^[a-zA-Z0-9_.]+(?=\\s*<-\\s*function)", txt, perl = TRUE))
  unique(m)
}

offending_exports <- function() {
  ns <- asNamespace("soilSIM")
  exports <- getNamespaceExports("soilSIM")
  fns <- exports[vapply(
    exports,
    function(x) is.function(get0(x, envir = ns, inherits = FALSE)),
    logical(1)
  )]
  bad <- fns[!grepl(approved_prefix_regex, fns)]
  setdiff(sort(bad), deprecated_alias_names())
}

test_that("all exported functions follow the approved verb/prefix lexicon", {
  skip("naming-convention enforcement is enabled in phase P6")
  expect_equal(offending_exports(), character(0))
})

test_that("public (non-internal) export count stays small", {
  skip("enabled once export demotion (phase P5) lands")
  # Target: <= 50 documented (non @keywords internal) exports.
  expect_lte(length(getNamespaceExports("soilSIM")), 60)
})
