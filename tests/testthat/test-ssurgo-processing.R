## Synthetic raw SSURGO-shaped tabular fixture (mimics mod01 download output:
## one row per horizon, component columns repeated/denormalized onto each
## horizon row, values held as character to exercise the string-parsing path).
make_raw_ssurgo_fixture <- function() {
  data.frame(
    cokey = c("1", "1", "2", "2", "3"),
    mukey = c("100", "100", "100", "100", "200"),
    compname = c("Alpha", "Alpha", "Beta", "Beta", "Gamma"),
    comppct_r = c(60, 60, 40, 40, 100),
    majcompflag = c("Yes", "Yes", "No", "No", "Yes"),
    taxclname = c("Fine, mixed", "Fine, mixed", NA, NA, "Coarse"),
    hzname = c("A", "Bw", "A", "Bt", "A"),
    hzdept_r = c(0, 20, 0, 15, 0),
    hzdepb_r = c(20, 45, 15, 60, 25),
    sandtotal_r = c("35-45", "30", "60", ">70", "TRACE"),
    claytotal_r = c("20", "25", "10", "8", "5"),
    silttotal_r = c("40", "45", "30", "22", "90"),
    dbovendry_r = c("1.3", "1.4", "1.1", "1.5", "1.2"),
    stringsAsFactors = FALSE
  )
}

test_that("identify_soil_property_columns_working() finds known SSURGO properties by their _r columns", {
  df <- make_raw_ssurgo_fixture()
  found <- identify_soil_property_columns_working(df)
  expect_setequal(found, c("sandtotal", "claytotal", "silttotal", "dbovendry"))

  df$unrelated_col_r <- 1
  expect_setequal(identify_soil_property_columns_working(df), c("sandtotal", "claytotal", "silttotal", "dbovendry"))
})

test_that("calculate_property_completeness_working() reports the non-NA fraction per property", {
  df <- make_raw_ssurgo_fixture()
  df$dbovendry_r <- suppressWarnings(as.numeric(df$dbovendry_r))
  df$dbovendry_r[2] <- NA
  completeness <- calculate_property_completeness_working(df)
  expect_equal(completeness$dbovendry, 4 / 5)
})

test_that("remove_invalid_horizons_working_compatible() drops bad-depth and missing-cokey rows only", {
  df <- data.frame(
    cokey = c("1", "2", "3", ""),
    hzdept_r = c(0, -5, 10, 0),
    hzdepb_r = c(20, 20, 5, 20),
    stringsAsFactors = FALSE
  )
  cleaned <- remove_invalid_horizons_working_compatible(df, max_depth = 250, verbose = FALSE)
  # row 1 valid; row 2 negative top depth; row 3 bottom <= top; row 4 empty cokey
  expect_equal(nrow(cleaned), 1)
  expect_equal(cleaned$cokey, "1")
})

test_that("remove_invalid_components_working() drops missing-cokey and out-of-range comppct_r rows only", {
  df <- data.frame(
    cokey = c("1", "", "3", "4"),
    comppct_r = c(60, 50, 150, NA),
    stringsAsFactors = FALSE
  )
  cleaned <- remove_invalid_components_working(df, verbose = FALSE)
  expect_setequal(cleaned$cokey, c("1", "4"))
})

test_that("calculate_component_stats_working() flags major components and taxonomic classification presence", {
  df <- data.frame(comppct_r = c(20, 10, NA), taxclname = c("Fine", NA, "Coarse"), stringsAsFactors = FALSE)
  result <- calculate_component_stats_working(df)
  expect_equal(result$is_major_component, c(TRUE, FALSE, FALSE))
  expect_equal(result$has_taxonomic_classification, c(TRUE, FALSE, TRUE))
})

test_that("calculate_derived_horizon_properties_working() fills thickness/midpoint/awc only when absent and computable", {
  df <- data.frame(hzdept_r = c(0, 20), hzdepb_r = c(20, 45), wthirdbar_r = c(30, 25), wfifteenbar_r = c(10, 8))
  result <- calculate_derived_horizon_properties_working(df)
  expect_equal(result$hzthk_r, c(20, 25))
  expect_equal(result$hz_midpoint, c(10, 32.5))
  expect_equal(result$awc_r, c(20, 17))

  # Existing hzthk_r is not overwritten.
  df2 <- data.frame(hzdept_r = 0, hzdepb_r = 20, hzthk_r = 999)
  expect_equal(calculate_derived_horizon_properties_working(df2)$hzthk_r, 999)
})

test_that("ensure_essential_columns_working() synthesizes cokey/hzname and coerces depth columns to numeric", {
  df <- data.frame(hzdept_r = c("0", "10"), hzdepb_r = c("10", "20"), stringsAsFactors = FALSE)
  result <- ensure_essential_columns_working(df)
  expect_true(all(c("cokey", "hzname") %in% names(result)))
  expect_true(is.numeric(result$hzdept_r))
  expect_true(is.numeric(result$hzdepb_r))
  expect_equal(result$hzname, c("Unknown", "Unknown"))
})

test_that("clean_property_data_ssurgo_compatible() parses character values via the consolidated canonical helpers", {
  df <- make_raw_ssurgo_fixture()
  result <- clean_property_data_ssurgo_compatible(df, "sandtotal", generate_report = FALSE)
  cleaned_vals <- result$data$sandtotal_r
  expect_true(is.numeric(cleaned_vals))
  # "35-45" -> range midpoint 40; "30" -> 30; "60" -> 60; ">70" -> 70*1.5=105
  # (then clamped out-of-range by apply_basic_range_limits below); "TRACE" -> 0.1
  expect_equal(cleaned_vals[1], 40)
  expect_equal(cleaned_vals[2], 30)
  expect_equal(cleaned_vals[5], 0.1)
})

test_that("clean_property_data_ssurgo_compatible() applies the union-merged range limits (mod02 partdensity/sand_* bounds now live in apply_basic_range_limits())", {
  df <- data.frame(cokey = "1", hzname = "A", hzdept_r = 0, hzdepb_r = 20,
                    partdensity_r = 5.0, sand_vc_r = 95, stringsAsFactors = FALSE)
  out_pd <- clean_property_data_ssurgo_compatible(df, "partdensity", generate_report = FALSE)$data
  # 5.0 is outside the merged (2.0, 3.2) bound -> clamped to NA by apply_basic_range_limits()
  expect_true(is.na(out_pd$partdensity_r))

  out_sand_vc <- clean_property_data_ssurgo_compatible(df, "sand_vc", generate_report = FALSE)$data
  # 95 is outside the merged (0, 80) bound -> clamped to NA
  expect_true(is.na(out_sand_vc$sand_vc_r))
})

test_that("process_ssurgo_data() runs end-to-end on the synthetic raw fixture and returns the documented shape", {
  raw <- make_raw_ssurgo_fixture()
  result <- process_ssurgo_data(raw, validate_results = TRUE, verbose = FALSE)

  expect_true(all(c("processed_data", "horizon_data", "component_data",
                     "processing_metadata", "validation_results", "quality_report") %in% names(result)))
  expect_true(nrow(result$processed_data) > 0)
  expect_true(nrow(result$component_data) > 0)
  expect_true(is.numeric(result$quality_report$overall_quality_score) ||
                is.na(result$quality_report$overall_quality_score))
  expect_true(result$quality_report$quality_grade %in%
                c("Excellent", "Good", "Acceptable", "Needs Improvement", "Poor", "Unknown"))
})

test_that("generate_processing_quality_report() renormalizes weights over available signals instead of zeroing out on missing ones", {
  full <- generate_processing_quality_report(
    original_data = data.frame(x = 1:10), processed_data = data.frame(x = 1:10),
    horizon_stats = list(property_completeness = list(a = 1)),
    component_stats = list(rows_retained = 1),
    validation_results = list(overall_quality = list(score = 1))
  )
  expect_equal(full$overall_quality_score, 1)

  partial <- generate_processing_quality_report(
    original_data = data.frame(x = 1:10), processed_data = data.frame(x = 1:10),
    horizon_stats = list(), component_stats = list(), validation_results = list()
  )
  # Only "retention" (=1) is available -> renormalized score is still 1, not
  # diluted toward 0 by the three missing signals.
  expect_equal(partial$overall_quality_score, 1)
})

test_that("hz_quant_prob_mukey() computes 05/50/95/PIW90 quantile columns per mukey/depth", {
  set.seed(1)
  hz_data <- data.frame(
    mukey = rep(c("1", "2"), each = 10),
    hzdept_r = 0, hzdepb_r = 20,
    db = c(stats::rnorm(10, 1.3, 0.1), stats::rnorm(10, 1.5, 0.1))
  )
  result <- hz_quant_prob_mukey(hz_data)

  expect_true(all(c("mukey", "top", "bottom", "Db_05", "Db_50", "Db_95", "Db_PIW90") %in% names(result)))
  expect_equal(nrow(result), 2)
  expect_true(all(result$Db_05 <= result$Db_50))
  expect_true(all(result$Db_50 <= result$Db_95))
  expect_true(all(result$Db_PIW90 >= 0))
})

test_that("hz_quant_prob_mukey() errors when no recognized property columns are present", {
  hz_data <- data.frame(mukey = "1", hzdept_r = 0, hzdepb_r = 20, not_a_property = 1)
  expect_error(hz_quant_prob_mukey(hz_data), "No valid soil properties")
})

test_that("hz_quant_prob_mukey() warns (not errors) on texture data without soiltexture installed", {
  testthat::skip_if(nzchar(system.file(package = "soiltexture")), "soiltexture is installed - texture-class path exercised separately")
  hz_data <- data.frame(
    mukey = "1", hzdept_r = 0, hzdepb_r = 20,
    sand_total = c(40, 42), silt_total = c(40, 38), clay_total = c(20, 20)
  )
  expect_warning(result <- hz_quant_prob_mukey(hz_data), "soiltexture package not available")
  expect_true("mukey" %in% names(result))
})

test_that("hz_quant_prob_mukey()'s single-pass quantile computation matches the original three-pass version", {
  # Regression test for the PERFORMANCE_IMPROVEMENT_PLAN.md Tier 3 hz_quant_prob_mukey() fix:
  # reimplements the ORIGINAL three-pass group_by()/summarize() + left_join() version's core
  # quantile/PIW90 computation here and checks the new single-pass version produces identical
  # numeric results, across multiple mukeys/depths/properties.
  set.seed(3)
  hz_data <- data.frame(
    mukey = rep(c("1", "2"), each = 20),
    hzdept_r = rep(c(0, 20), 20),
    hzdepb_r = rep(c(20, 40), 20),
    db = c(stats::rnorm(20, 1.3, 0.1), stats::rnorm(20, 1.5, 0.15)),
    ph = c(stats::rnorm(20, 6.0, 0.3), stats::rnorm(20, 6.5, 0.2))
  )

  old_three_pass <- function(hz_data) {
    q <- c(0.05, 0.5, 0.95)
    data_out <- hz_data |>
      dplyr::select(mukey, top = hzdept_r, bottom = hzdepb_r, Db = db, ph = ph)
    prop_names <- c("Db", "ph")

    data_stats05 <- data_out |> dplyr::group_by(mukey, top) |>
      dplyr::summarize(dplyr::across(dplyr::all_of(prop_names), ~stats::quantile(.x, probs = q[1], na.rm = TRUE)), .groups = "drop") |>
      as.data.frame()
    data_stats50 <- data_out |> dplyr::group_by(mukey, top) |>
      dplyr::summarize(dplyr::across(dplyr::all_of(prop_names), ~stats::quantile(.x, probs = q[2], na.rm = TRUE)), .groups = "drop") |>
      as.data.frame()
    data_stats95 <- data_out |> dplyr::group_by(mukey, top) |>
      dplyr::summarize(dplyr::across(dplyr::all_of(prop_names), ~stats::quantile(.x, probs = q[3], na.rm = TRUE)), .groups = "drop") |>
      as.data.frame()

    data_statsPIW90 <- data_stats95 |> dplyr::select(-mukey, -top) - data_stats05 |> dplyr::select(-mukey, -top)
    data_statsPIW90 <- data_statsPIW90 |>
      dplyr::mutate(mukey = data_stats05$mukey, .before = 1) |>
      dplyr::mutate(top = data_stats05$top, .after = mukey)

    names(data_stats05)[-(1:2)] <- paste0(names(data_stats05)[-(1:2)], "_05")
    names(data_stats50)[-(1:2)] <- paste0(names(data_stats50)[-(1:2)], "_50")
    names(data_stats95)[-(1:2)] <- paste0(names(data_stats95)[-(1:2)], "_95")
    names(data_statsPIW90)[-(1:2)] <- paste0(names(data_statsPIW90)[-(1:2)], "_PIW90")

    data_stats05 |>
      dplyr::left_join(data_stats50, by = c("mukey", "top")) |>
      dplyr::left_join(data_stats95, by = c("mukey", "top")) |>
      dplyr::left_join(data_statsPIW90, by = c("mukey", "top"))
  }

  expected <- old_three_pass(hz_data)
  actual <- hz_quant_prob_mukey(hz_data)
  common_cols <- c("mukey", "top", "Db_05", "Db_50", "Db_95", "Db_PIW90", "ph_05", "ph_50", "ph_95", "ph_PIW90")
  expect_equal(
    actual[order(actual$mukey, actual$top), common_cols],
    expected[order(expected$mukey, expected$top), common_cols],
    tolerance = 1e-12
  )
})
