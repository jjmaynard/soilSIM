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

test_that("fetch_ssurgo_all_components_working() requires the live SDA service", {
  testthat::skip_if_offline()
  testthat::skip("Live NRCS Soil Data Access queries are not exercised in automated tests.")
})

test_that("synthesize_component_horizons_from_siblings() averages numeric properties per genhz group across siblings", {
  fx <- make_component_recovery_fixture()
  target <- fx$all_components[fx$all_components$cokey == "20", , drop = FALSE]
  siblings <- fx$ssurgo_data[fx$ssurgo_data$compname == "compA", , drop = FALSE]

  result <- synthesize_component_horizons_from_siblings(target, siblings)

  expect_equal(nrow(result), 2)  # one row per genhz group (A, Bt)
  expect_setequal(classify_genhz(result$hzname), c("A", "B"))
  expect_true(all(result$infill_method == "component_aoi_average"))
  expect_true(all(result$component_synthesized))
  # identity columns come from the missing component, never averaged from siblings
  expect_true(all(result$mukey == target$mukey))
  expect_true(all(result$cokey == target$cokey))
  expect_true(all(result$compname == target$compname))
  expect_true(all(result$comppct_r == target$comppct_r))

  a_row <- result[classify_genhz(result$hzname) == "A", ]
  expect_equal(a_row$sandtotal_r, mean(c(40, 44)))
  expect_equal(a_row$claytotal_r, mean(c(15, 19)))
  bt_row <- result[classify_genhz(result$hzname) == "B", ]
  expect_equal(bt_row$sandtotal_r, mean(c(30, 34)))
  expect_equal(bt_row$claytotal_r, mean(c(25, 29)))

  # chkey is freshly generated and unique, never borrowed from a sibling
  expect_true(all(grepl("^synth_20_", result$chkey)))
})

test_that("synthesize_component_horizons_from_siblings() with a single sibling is a copy (mean of one)", {
  fx <- make_component_recovery_fixture()
  target <- fx$all_components[fx$all_components$cokey == "20", , drop = FALSE]
  one_sibling <- fx$ssurgo_data[fx$ssurgo_data$cokey == "10", , drop = FALSE]  # compA, sibling1 only

  result <- synthesize_component_horizons_from_siblings(target, one_sibling)

  expect_equal(nrow(result), 2)
  a_row <- result[classify_genhz(result$hzname) == "A", ]
  expect_equal(a_row$sandtotal_r, 40)
  expect_equal(a_row$claytotal_r, 15)
})

test_that("synthesize_component_horizons_from_siblings() still produces a genhz group present in only some siblings", {
  fx <- make_component_recovery_fixture()
  target <- fx$all_components[fx$all_components$cokey == "20", , drop = FALSE]
  siblings <- fx$ssurgo_data[fx$ssurgo_data$compname == "compA", , drop = FALSE]
  # Give sibling2 a Cr horizon that sibling1 doesn't have, at the record level (append, don't replace).
  extra_cr <- siblings[1, , drop = FALSE]
  extra_cr$cokey <- "11"
  extra_cr$hzname <- "Cr"
  extra_cr$hzdept_r <- 55
  extra_cr$hzdepb_r <- 100
  extra_cr$sandtotal_r <- 20
  extra_cr$claytotal_r <- 10
  siblings <- rbind(siblings, extra_cr)

  result <- synthesize_component_horizons_from_siblings(target, siblings)

  expect_setequal(classify_genhz(result$hzname), c("A", "B", "Cr"))
  cr_row <- result[classify_genhz(result$hzname) == "Cr", ]
  expect_equal(nrow(cr_row), 1)
  expect_equal(cr_row$sandtotal_r, 20)  # averaged over exactly the 1 sibling that had it
})

test_that("synthesize_component_horizons_from_siblings() returns 0 rows when no sibling horizon is classifiable", {
  fx <- make_component_recovery_fixture()
  target <- fx$all_components[fx$all_components$cokey == "20", , drop = FALSE]
  siblings <- fx$ssurgo_data[fx$ssurgo_data$compname == "compA", , drop = FALSE]
  siblings$hzname <- "H1"  # unparseable by classify_genhz()

  result <- synthesize_component_horizons_from_siblings(target, siblings)
  expect_equal(nrow(result), 0)
})

test_that("recover_missing_horizon_components() recovers a component via an AOI-wide sibling search (different mukey)", {
  fx <- make_component_recovery_fixture()
  result <- recover_missing_horizon_components(fx$ssurgo_data, fx$all_components, verbose = FALSE)

  expect_equal(nrow(result$components_missing_horizons), 0)
  expect_equal(nrow(result$components_recovered), 1)
  expect_equal(result$components_recovered$cokey, "20")
  expect_equal(result$components_recovered$n_sibling_cokeys_used, 2)

  recovered_rows <- result$ssurgo_data[result$ssurgo_data$cokey == "20", ]
  expect_equal(nrow(recovered_rows), 2)
  expect_true(all(recovered_rows$component_synthesized))

  # Provenance backfill: no row anywhere (real or synthesized) has NA in either column
  expect_false(any(is.na(result$ssurgo_data$infill_method)))
  expect_false(any(is.na(result$ssurgo_data$component_synthesized)))
  # Real rows are untouched (not flagged as synthesized)
  real_rows <- result$ssurgo_data[result$ssurgo_data$cokey != "20", ]
  expect_true(all(!real_rows$component_synthesized))
  expect_true(all(real_rows$infill_method == ""))
})

test_that("recover_missing_horizon_components() leaves a component unresolved when no AOI sibling shares its compname", {
  fx <- make_component_recovery_fixture()
  # A component with a compname that appears nowhere else in ssurgo_data.
  lonely <- data.frame(mukey = "1", cokey = "30", compname = "compZ",
                        comppct_l = 5, comppct_r = 8, comppct_h = 11, stringsAsFactors = FALSE)
  all_components <- rbind(fx$all_components, lonely)

  result <- recover_missing_horizon_components(fx$ssurgo_data, all_components, verbose = FALSE)

  expect_equal(nrow(result$ssurgo_data), nrow(fx$ssurgo_data) + 2)  # only cokey 20 (compA) recovered
  expect_true("30" %in% result$components_missing_horizons$cokey)
  expect_false("30" %in% result$components_recovered$cokey)
})

test_that("recover_missing_horizon_components() is a no-op when every real component already has horizon data", {
  fx <- make_component_recovery_fixture()
  all_present <- fx$all_components[fx$all_components$cokey != "20", ]  # drop the missing one

  result <- recover_missing_horizon_components(fx$ssurgo_data, all_present, verbose = FALSE)

  expect_identical(result$ssurgo_data, fx$ssurgo_data)
  expect_equal(nrow(result$components_missing_horizons), 0)
  expect_equal(nrow(result$components_recovered), 0)
})

test_that("recover_missing_horizon_components() end-to-end: a recovered component gets a real sim_comppct via sim_component_comp()", {
  fx <- make_component_recovery_fixture()
  result <- recover_missing_horizon_components(fx$ssurgo_data, fx$all_components, verbose = FALSE)

  hz_data <- result$ssurgo_data
  hz_data$genhz <- classify_genhz(hz_data$hzname)
  expect_false(any(is.na(hz_data$genhz)))

  comp_draws <- sim_component_comp(hz_data, n_simulations = 50)
  recovered_comp <- comp_draws[comp_draws$cokey == "20", ]
  expect_equal(nrow(recovered_comp), 1)
  expect_false(is.na(recovered_comp$sim_comppct))
  expect_true(recovered_comp$sim_comppct > 0)
})

test_that("create_download_metadata() reports component-recovery counts, and NA when not supplied", {
  base_args <- list(
    start_time = Sys.time() - 1, end_time = Sys.time(), aoi_wkt = sample_wkt,
    properties = c("sandtotal"), include_restrictions = TRUE,
    mukey_count = 1, data_rows = 10, unique_cokeys = 2,
    spatial_result = list(aoi = NULL), validation_results = NULL
  )

  fx <- make_component_recovery_fixture()
  with_counts <- do.call(create_download_metadata, c(base_args, list(
    components_missing_horizons = fx$all_components[0, ],
    components_recovered = data.frame(mukey = "2", cokey = "20", compname = "compA",
                                       n_sibling_cokeys_used = 2, stringsAsFactors = FALSE)
  )))
  expect_equal(with_counts$components_missing_horizons_count, 0)
  expect_equal(with_counts$components_recovered_count, 1)

  without_counts <- do.call(create_download_metadata, base_args)
  expect_true(is.na(without_counts$components_missing_horizons_count))
  expect_true(is.na(without_counts$components_recovered_count))
})
