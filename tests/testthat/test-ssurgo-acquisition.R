## Only pure-computation helpers are exercised directly here; anything that
## hits the live SDA web service (download_ssurgo_tabular(),
## process_aoi_and_get_mukeys_working(), execute_ssurgo_query_working()) is
## skipped rather than mocked - SDA's response shape is too fragile to
## fixture-maintain for the value it'd add.
sample_wkt <- "POLYGON((-120.5 38.5, -120.4 38.5, -120.4 38.6, -120.5 38.6, -120.5 38.5))"

test_that("create_ssurgo_property_lookup_working() returns the 9-property low/rep/high lookup table", {
  lookup <- create_ssurgo_property_lookup_working()
  expect_setequal(lookup$Property, c("sandtotal", "claytotal", "silttotal", "dbovendry",
                                      "ph1to1h2o", "cec7", "om", "wthirdbar", "wfifteenbar"))
  expect_true(all(c("SSURGO_Label_Low", "SSURGO_Label_Rep", "SSURGO_Label_High") %in% names(lookup)))
  expect_equal(lookup$SSURGO_Label_Rep[lookup$Property == "sandtotal"], "sandtotal_r")
})

test_that("REGRESSION: get_predefined_properties('ssurgo') now resolves the real lookup instead of character(0)", {
  # Before mod01's migration into this package, create_ssurgo_property_lookup_working()
  # did not exist, so utils.R's tryCatch() silently degraded to character(0).
  props <- get_predefined_properties("ssurgo")
  expect_length(props, 9)
  expect_true("dbovendry" %in% props)
})

test_that("aggregate_rock_fragment_volume_working() sums RFV across fragment-size classes within a horizon", {
  df <- data.frame(
    chkey = c("h1", "h1", "h2"),
    fragsize_r = c(5, 20, 10),
    rfv_l = c(1, 2, 3), rfv_r = c(2, 3, 4), rfv_h = c(3, 4, 5),
    stringsAsFactors = FALSE
  )
  result <- aggregate_rock_fragment_volume_working(df, verbose = FALSE)
  # The function replaces the per-fragsize RFV columns with the horizon-level
  # aggregate but does NOT collapse fragsize_r duplicate rows away (it
  # left_join()s the one-row-per-chkey aggregate back onto the original,
  # still-fragsize-duplicated rows) - so h1 keeps its original 2 rows, both
  # now carrying the same summed rfv_r (2 + 3 = 5).
  h1 <- result[result$chkey == "h1", ]
  expect_equal(nrow(h1), 2)
  expect_equal(h1$rfv_r, c(5, 5))

  h2 <- result[result$chkey == "h2", ]
  expect_equal(h2$rfv_r, 4)
})

test_that("aggregate_rock_fragment_volume_working() is a no-op when RFV or fragsize columns are absent", {
  df <- data.frame(chkey = "h1", hzname = "A", stringsAsFactors = FALSE)
  expect_identical(aggregate_rock_fragment_volume_working(df, verbose = FALSE), df)
})

test_that("add_restriction_indicators_working() flags Cr/Cd/Cx horizon names and O/R master designations", {
  df <- data.frame(
    hzname = c("Bt", "2Cr", "Oi", "R"),
    desgnmaster = c("B", NA, "O", "R"),
    hzdept_r = c(0, 60, 0, 100),
    dbovendry_r = c(1.3, 1.4, 1.2, 1.5),
    reskind = c(NA, "lithic bedrock", NA, "paralithic bedrock"),
    resdept_l = c(NA, 55, NA, 95), resdept_r = c(NA, 60, NA, 100),
    texcl = NA_character_, lieutex = NA_character_,
    stringsAsFactors = FALSE
  )
  result <- add_restriction_indicators_working(df, verbose = FALSE)
  expect_true(result$horizon_suggests_restriction[2])   # "2Cr" matches Cr pattern
  expect_true(result$organic_horizon[3])                # desgnmaster == "O"
  expect_true(result$horizon_suggests_restriction[4])   # desgnmaster == "R"
  expect_false(result$organic_horizon[1])
  expect_true(is.logical(result$unsuitable_horizon))
})

test_that("generate_ssurgo_cache_key() is deterministic and order-independent in properties, but sensitive to AOI/restriction flag", {
  key1 <- generate_ssurgo_cache_key(sample_wkt, c("sandtotal", "claytotal"), TRUE)
  key2 <- generate_ssurgo_cache_key(sample_wkt, c("claytotal", "sandtotal"), TRUE)
  expect_identical(key1, key2)

  key_no_restrictions <- generate_ssurgo_cache_key(sample_wkt, c("sandtotal", "claytotal"), FALSE)
  expect_false(identical(key1, key_no_restrictions))

  key_other_aoi <- generate_ssurgo_cache_key("POLYGON((0 0,1 0,1 1,0 1,0 0))", c("sandtotal", "claytotal"), TRUE)
  expect_false(identical(key1, key_other_aoi))

  expect_match(key1, "^ssurgo_[A-Za-z0-9_-]+$")
})

test_that("validate_download_inputs_ssurgo() passes for a well-formed AOI/properties and fails when a required parameter is missing", {
  # suppressWarnings(): local GDAL/PROJ install version-mismatch warnings from
  # terra (environment-specific, unrelated to this function's logic).
  valid_result <- suppressWarnings(validate_download_inputs_ssurgo(
    aoi_wkt = sample_wkt,
    properties = c("sandtotal", "claytotal"),
    include_restrictions = TRUE
  ))
  expect_true(valid_result$valid)
  expect_equal(length(valid_result$errors), 0)

  # properties is a required parameter (per param_specs) - character(0) triggers
  # the "Required parameter missing" branch and an early return, deterministically,
  # without depending on how a malformed WKT string happens to be handled
  # downstream (validate_wkt_geometry()'s own error path - a pre-existing utils.R
  # quirk out of scope for this migration - logs a WARN rather than an error and
  # does not flip `valid` to FALSE; not exercised here).
  invalid_result <- validate_download_inputs_ssurgo(
    aoi_wkt = sample_wkt,
    properties = character(0),
    include_restrictions = TRUE
  )
  expect_false(invalid_result$valid)
  expect_true(length(invalid_result$errors) > 0)
})

test_that("check_ssurgo_cache()/cache_ssurgo_data() round-trip through a tempdir", {
  cache_dir <- file.path(tempdir(), paste0("ssurgo_cache_test_", as.integer(Sys.time())))
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)
  data <- data.frame(cokey = c("1", "2"), hzdept_r = c(0, 10), hzdepb_r = c(10, 20))
  mu <- NULL
  properties <- c("sandtotal", "claytotal")

  miss <- check_ssurgo_cache(sample_wkt, properties, TRUE, cache_dir, verbose = FALSE)
  expect_null(miss)

  cache_result <- cache_ssurgo_data(data, mu, sample_wkt, properties, TRUE, cache_dir, verbose = FALSE)
  expect_true(cache_result$cached)
  expect_true(file.exists(cache_result$cache_file))

  hit <- check_ssurgo_cache(sample_wkt, properties, TRUE, cache_dir, verbose = FALSE)
  expect_true(hit$cache_hit)
  expect_equal(hit$data, data)
})

test_that("download_ssurgo_tabular()/process_aoi_and_get_mukeys_working()/execute_ssurgo_query_working() require the live SDA service", {
  testthat::skip_if_offline()
  testthat::skip("Live NRCS Soil Data Access queries are not exercised in automated tests - see B.6 in the migration plan.")
})
