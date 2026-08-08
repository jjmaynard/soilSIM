test_that("parse_single_string_advanced() parses ranges, comparisons, and qualitative terms", {
  expect_equal(parse_single_string_advanced("15-20")$value, 17.5)
  expect_equal(parse_single_string_advanced("15-20")$type, "range")
  # NOTE: "TO" in the range regex's character class matches single T/O
  # characters, not the literal word "TO" - so "10 TO 15" does NOT hit the
  # range branch and instead falls through to a plain leading-number
  # extraction. This is pre-existing behavior carried over unchanged from
  # modules/mod03_data_infilling.R, not something introduced by this
  # migration - locking in the actual (if surprising) current behavior here.
  expect_equal(parse_single_string_advanced("10 TO 15")$value, 10)
  expect_equal(parse_single_string_advanced("10 TO 15")$type, "numeric")
  # Non-ASCII dash characters (en dash, em dash) - regex now expressed via
  # \uXXXX escapes in the source, must still match at runtime.
  expect_equal(parse_single_string_advanced("5–20")$value, 12.5)  # en dash
  expect_equal(parse_single_string_advanced("5—20")$value, 12.5)  # em dash
  expect_equal(parse_single_string_advanced("<5")$type, "upper_bound")
  expect_equal(parse_single_string_advanced(">50")$type, "lower_bound")
  # <= / >= unicode comparison operators
  expect_equal(parse_single_string_advanced("≤10")$type, "upper_bound")  # <=10
  expect_equal(parse_single_string_advanced("≥25")$type, "lower_bound")  # >=25
  expect_equal(parse_single_string_advanced("TRACE")$type, "qualitative")
  expect_equal(parse_single_string_advanced("TRACE")$value, 0.1)
})

test_that("calculate_saxton_rawls_single() produces physically plausible, checkable output", {
  result <- calculate_saxton_rawls_single(sand_pct = 40, clay_pct = 20, silt_pct = 40, bulk_density = 1.3)
  expect_true(result$field_capacity > result$wilting_point)
  expect_true(result$available_water_capacity > 0)
  expect_equal(result$available_water_capacity, result$field_capacity - result$wilting_point, tolerance = 0.5)
  expect_true(result$field_capacity_l < result$field_capacity)
  expect_true(result$field_capacity_h > result$field_capacity)

  # Texture renormalization: percentages summing far from 100 still produce
  # finite, plausible output (rescaled internally).
  result_bad_texture <- calculate_saxton_rawls_single(sand_pct = 50, clay_pct = 50, silt_pct = 50, bulk_density = 1.3)
  expect_true(all(is.finite(unlist(result_bad_texture))))
})

test_that("impute_rfv_values() texture-informed default and valid-value spread behave as documented", {
  missing_rfv <- impute_rfv_values(list(claytotal_r = 20, sandtotal_r = 40, silttotal_r = 40, rfv_r = NA))
  expect_equal(missing_rfv$rfv_r, 0.5)
  expect_true(missing_rfv$rfv_l < missing_rfv$rfv_r && missing_rfv$rfv_r < missing_rfv$rfv_h)

  valid_rfv <- impute_rfv_values(list(claytotal_r = 20, sandtotal_r = 40, silttotal_r = 40, rfv_r = 20))
  expect_equal(valid_rfv$rfv_r, 20)
  expect_equal(valid_rfv$rfv_l, 14, tolerance = 1e-6)
  expect_equal(valid_rfv$rfv_h, 26, tolerance = 1e-6)

  zero_rfv <- impute_rfv_values(list(claytotal_r = 20, sandtotal_r = 40, silttotal_r = 40, rfv_r = 0))
  expect_equal(zero_rfv$rfv_r, 0.1)
})

test_that("cross_component_property_interpolation() fills from another component at a similar depth", {
  group <- data.frame(
    cokey = c("1", "2"),
    hzname = c("A", "A"),
    hzdept_r = c(0, 0),
    hzdepb_r = c(20, 15),
    om_r = c(2, NA),
    unsuitable_horizon = c(FALSE, FALSE)
  )
  res <- cross_component_property_interpolation(group, "om_r")
  # Row 2 (cokey "2") has a depth midpoint of 7.5cm, row 1's midpoint is 10cm -
  # within the 15cm tolerance, so row 2 borrows row 1's om_r.
  expect_equal(res$om_r[2], 2)
  expect_true(grepl("cross_comp", res$infill_method[2]))
  expect_equal(res$om_r[1], 2)  # row 1 (the source) is untouched
})

test_that("cross_component_property_interpolation() leaves a value NA when no source is within depth tolerance", {
  group <- data.frame(
    cokey = c("1", "2"),
    hzname = c("A", "A"),
    hzdept_r = c(0, 100),
    hzdepb_r = c(20, 120),
    om_r = c(2, NA),
    unsuitable_horizon = c(FALSE, FALSE)
  )
  res <- cross_component_property_interpolation(group, "om_r")
  expect_true(is.na(res$om_r[2]))  # depth midpoints are 10cm vs 110cm - far outside the 15cm tolerance
})

test_that("related_property_estimation() applies the texture sum-to-100 constraint", {
  group <- make_horizon_row(sandtotal_r = 40, silttotal_r = 35)
  group$claytotal_r <- NA_real_
  group$unsuitable_horizon <- FALSE
  config <- get_default_property_config("claytotal")
  res <- related_property_estimation(group, "claytotal", config)
  expect_equal(res$claytotal_r, 25)  # 100 - 40 - 35
})

test_that("related_property_estimation() applies clay-based water retention formulas", {
  group <- make_horizon_row(claytotal_r = 30)
  group$wthirdbar_r <- NA_real_
  group$unsuitable_horizon <- FALSE
  config <- get_default_property_config("wthirdbar")
  res <- related_property_estimation(group, "wthirdbar", config)
  expect_equal(res$wthirdbar_r, 0.3 * 30 + 10)
})

test_that("related_property_estimation() applies clay/OM-based CEC estimation", {
  # "om" isn't one of make_horizon_row()'s default properties, so its _r
  # column must be added explicitly rather than passed as an override.
  group <- make_horizon_row(claytotal_r = 20)
  group$om_r <- 3
  group$cec7_r <- NA_real_
  group$unsuitable_horizon <- FALSE
  config <- get_default_property_config("cec7")
  res <- related_property_estimation(group, "cec7", config)
  expect_equal(res$cec7_r, max(2, 20 * 0.5 + 3 * 20))
})

test_that("related_property_estimation() applies horizon/OM-adjusted pH estimation", {
  group <- make_horizon_row()
  group$om_r <- 6  # "om" isn't a default property - add its _r column explicitly
  group$ph1to1h2o_r <- NA_real_
  group$unsuitable_horizon <- FALSE
  # Like dbovendry (see below), get_default_property_config("ph1to1h2o") does
  # not set related_properties (pre-existing legacy config characteristic),
  # so supply it explicitly to exercise the branch directly.
  config <- get_default_property_config("ph1to1h2o")
  config$related_properties <- c("om")
  res <- related_property_estimation(group, "ph1to1h2o", config)
  expect_equal(res$ph1to1h2o_r, 6.2 - 0.3 - 0.4)
})

test_that("related_property_estimation() is a no-op for ph1to1h2o's actual default config (related_properties unset)", {
  group <- make_horizon_row()
  group$om_r <- 6
  group$ph1to1h2o_r <- NA_real_
  group$unsuitable_horizon <- FALSE
  config <- get_default_property_config("ph1to1h2o")
  res <- related_property_estimation(group, "ph1to1h2o", config)
  expect_true(is.na(res$ph1to1h2o_r))
})

test_that("related_property_estimation() applies depth/horizon/clay-adjusted organic matter estimation", {
  group <- make_horizon_row(hzname = "A", hzdept_r = 5, claytotal_r = 40)
  group$om_r <- NA_real_
  group$unsuitable_horizon <- FALSE
  config <- get_default_property_config("om")
  res <- related_property_estimation(group, "om", config)
  expect_equal(res$om_r, 3.5 * 1.5 * 1.3)
})

test_that("related_property_estimation() applies texture-adjusted bulk density estimation", {
  group <- make_horizon_row(claytotal_r = 30, sandtotal_r = 40)
  group$dbovendry_r <- NA_real_
  group$unsuitable_horizon <- FALSE
  # get_default_property_config("dbovendry") does not set related_properties
  # (a pre-existing characteristic of the legacy reference config, carried
  # over unchanged - not something introduced by this port), so
  # related_property_estimation() bails out before reaching the
  # bulk_density branch when called with the *default* config. Supplying
  # related_properties explicitly here exercises that branch directly.
  config <- get_default_property_config("dbovendry")
  config$related_properties <- c("claytotal", "sandtotal")
  res <- related_property_estimation(group, "dbovendry", config)
  expect_equal(res$dbovendry_r, 1.4 - (30 - 20) * 0.01 + (40 - 50) * 0.005)
})

test_that("related_property_estimation() is a no-op for dbovendry's actual default config (related_properties unset)", {
  group <- make_horizon_row(claytotal_r = 30, sandtotal_r = 40)
  group$dbovendry_r <- NA_real_
  group$unsuitable_horizon <- FALSE
  config <- get_default_property_config("dbovendry")
  res <- related_property_estimation(group, "dbovendry", config)
  expect_true(is.na(res$dbovendry_r))
})

test_that("infill_soil_property() fills missing _l/_h from a complete horizon, and recovers a value via cross-component interpolation", {
  df <- data.frame(
    cokey = c("1", "1", "2"),
    hzname = c("A", "B", "A"),
    hzdept_r = c(0, 20, 0),
    hzdepb_r = c(20, 40, 15),
    om_l = c(1, NA, 2), om_r = c(2, 3, NA), om_h = c(4, NA, 6)
  )
  res <- infill_soil_property(df, "om")
  expect_true(is.finite(res$om_h[2]))  # row 2's om_h recovered from its own om_r + a learned/contextual spread
  expect_true(is.finite(res$om_l[2]))
  # Row 3 (cokey "2") has a depth midpoint of 7.5cm; row 1 (cokey "1") has a
  # midpoint of 10cm - within cross-component interpolation's 15cm tolerance,
  # so row 3's om_r is now recovered from row 1's om_r (=2) rather than
  # staying NA. This is the intended behavior once cross-component
  # interpolation (Strategy 4) is implemented for real - previously, with
  # that strategy stubbed out, row 3 had no way to borrow from another
  # component and stayed NA.
  expect_equal(res$om_r[3], 2)
})

test_that("get_default_property_config()/create_custom_property_config() return usable configs", {
  db_config <- get_default_property_config("dbovendry")
  expect_equal(db_config$units, "g/cm^3")  # non-ASCII unit symbol replaced with ASCII-safe text
  expect_true(is.numeric(db_config$typical_range))

  custom <- create_custom_property_config("my_property", property_type = "generic", typical_range = c(0, 10))
  expect_equal(custom$typical_range, c(0, 10))
})

test_that("validate_property_config() accepts a well-formed config and rejects malformed ones", {
  good <- list(type = "generic", units = "%", fallback_range = 5)
  expect_true(validate_property_config(good, "test_prop"))

  expect_error(validate_property_config(list(type = "generic"), "test_prop"), "missing required fields")
  expect_error(validate_property_config(list(type = "generic", units = "%", fallback_range = -1), "test_prop"),
               "must be a positive number")
  expect_error(validate_property_config(list(type = "generic", units = "%", fallback_range = 1,
                                             typical_range = c(10, 5)), "test_prop"),
               "min < max")
  expect_error(validate_property_config("not a list", "test_prop"), "must be a list")
})

test_that("get_rfv_range_category() categorizes RFV values at documented boundaries", {
  expect_equal(get_rfv_range_category(NA), "none")
  expect_equal(get_rfv_range_category(0), "none")
  expect_equal(get_rfv_range_category(5), "low")
  expect_equal(get_rfv_range_category(5.001), "moderate")
  expect_equal(get_rfv_range_category(15), "moderate")
  expect_equal(get_rfv_range_category(15.001), "high")
  expect_equal(get_rfv_range_category(35), "high")
  expect_equal(get_rfv_range_category(35.001), "very_high")
  expect_equal(get_rfv_range_category(60), "very_high")
  expect_equal(get_rfv_range_category(61), "extreme")
})

test_that("apply_property_constraints() clamps per property type and typical_range", {
  texture_config <- list(type = "texture", units = "%", fallback_range = 5, typical_range = c(0, 100))
  expect_equal(apply_property_constraints(c(-10, 50, 150), texture_config), c(0, 50, 100))

  ph_config <- list(type = "ph", units = "pH units", fallback_range = 0.5)
  expect_equal(apply_property_constraints(c(-1, 7, 20), ph_config), c(0, 7, 14))

  rfv_config <- list(type = "rock_fragments", units = "%", fallback_range = 5)
  expect_equal(apply_property_constraints(c(-5, 50, 200), rfv_config), c(0.01, 50, 95))

  generic_config <- list(type = "generic", units = "unknown", fallback_range = 5)
  expect_equal(apply_property_constraints(c(-3, 10), generic_config), c(0, 10))

  expect_true(all(is.na(apply_property_constraints(c(NA, NA), texture_config))))
})

test_that("create_validation_config()/add_range_rule()/apply_validation_rules() round-trip correctly", {
  config <- create_validation_config()
  expect_equal(config$range_rules, list())

  config <- add_range_rule(config, "sandtotal", min_val = 0, max_val = 100)
  result <- apply_validation_rules(c(-5, 50, 150), "sandtotal", config)
  expect_equal(result$violations, c(TRUE, FALSE, TRUE))
  expect_equal(result$rules_applied, 1)

  # A rule for a different property doesn't apply
  result_other <- apply_validation_rules(c(-5, 50, 150), "claytotal", config)
  expect_equal(result_other$violations, c(FALSE, FALSE, FALSE))
})

test_that("add_relationship_rule() records the rule but apply_validation_rules() does not yet enforce it (matches legacy behavior)", {
  config <- create_validation_config()
  config <- add_relationship_rule(config, c("sandtotal", "claytotal", "silttotal"),
                                  relationship_type = "sum", expected_sum = 100)
  expect_equal(length(config$relationship_rules), 1)
  # No range_rules were added, so nothing is flagged even for a clearly
  # invalid value - relationship_rules are recorded but not consumed, as in
  # the original reference implementation.
  result <- apply_validation_rules(c(1000), "sandtotal", config)
  expect_equal(result$violations, FALSE)
  expect_equal(result$rules_applied, 0)
})

test_that("summarize_unsuitable_horizons() summarizes excluded horizon types", {
  df <- data.frame(
    hzname = c("A", "O", "R", "Bt"),
    unsuitable_horizon = c(FALSE, TRUE, TRUE, FALSE)
  )
  result <- summarize_unsuitable_horizons(df)
  expect_equal(result$n_unsuitable, 2)
  expect_setequal(result$horizon_types, c("O", "R"))

  expect_error(summarize_unsuitable_horizons(data.frame(x = 1)), "unsuitable_horizon")
})

test_that("infill_property_range_values() classifies depth zones correctly (regression for an apply()-coercion bug)", {
  # Regression test for the PERFORMANCE_IMPROVEMENT_PLAN.md Tier 2 infill_property_range_values()
  # fix. The OLD apply(df, 1, ...) implementation coerced the entire row (all of df's columns) to
  # character via as.matrix() before calling calculate_property_lower_bound()/upper_bound(), which
  # read hzdepb_r WITHOUT as.numeric() for its depth-zone classification (if (depth <= 30) ...) -
  # so with hzdepb_r coerced to character, "5" <= "30" was a STRING comparison (FALSE, since "5" >
  # "30" lexicographically), silently misclassifying a 5cm-deep horizon as "depth_deep" instead of
  # "depth_surface". The new column-subsetted, type-preserving implementation fixes this. Distinct
  # hznames per row (never repeating 3x) keep the horizon-specific learned range from ever
  # qualifying, isolating the depth-zone branch as the one that determines the spread used.
  surface_rows <- data.frame(
    hzname = c("Ax1", "Ax2", "Ax3"), hzdepb_r = c(10, 20, 25),
    testprop_r = 10, testprop_l = 8, testprop_h = 12  # spread = 2 both sides
  )
  deep_rows <- data.frame(
    hzname = c("Cx1", "Cx2", "Cx3"), hzdepb_r = c(120, 150, 180),
    testprop_r = 10, testprop_l = 0, testprop_h = 20  # spread = 10 both sides
  )
  target_row <- data.frame(
    hzname = "Ax_target", hzdepb_r = 5,  # shallow - must classify as depth_surface, not depth_deep
    testprop_r = 10, testprop_l = NA_real_, testprop_h = NA_real_
  )
  df <- rbind(surface_rows, deep_rows, target_row)

  property_config <- list(type = "other", fallback_range = 99, typical_range = NULL)
  result <- infill_property_range_values(df, "testprop", property_config)

  target_idx <- nrow(result)
  # Correct (depth_surface) classification -> spread 2, not the depth_deep spread of 10.
  expect_equal(result$testprop_l[target_idx], 8)
  expect_equal(result$testprop_h[target_idx], 12)
})
