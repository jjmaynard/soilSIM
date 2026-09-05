test_that("classify_genhz() correctly distinguishes Cr from C1 and handles leading digit prefixes", {
  result <- classify_genhz(c("Bt2", "2Bt3", "Cr", "C1", "2C", "Oi", "R", "AB", "BC", NA))
  expect_equal(result, c("B", "B", "Cr", "C", "C", "O", "R", "A", "B", NA))
})

test_that("classify_genhz() returns NA for unrecognized hzname text, not an error", {
  expect_equal(classify_genhz("Xyz123"), NA_character_)
  expect_equal(classify_genhz(""), NA_character_)
  expect_no_error(classify_genhz(character(0)))
  expect_length(classify_genhz(character(0)), 0)
})

test_that("classify_genhz() is vectorized and case-insensitive", {
  expect_equal(classify_genhz(c("a", "bt", "cr")), c("A", "B", "Cr"))
})

test_that("build_kssl_fallback_matrix() returns a PD block-diagonal matrix with identity for unmapped properties", {
  result <- build_kssl_fallback_matrix(c("dbovendry", "ph1to1h2o", "unmapped_prop"), genhz = "A")
  expect_equal(result["dbovendry", "unmapped_prop"], 0)
  expect_equal(result["unmapped_prop", "unmapped_prop"], 1)
  expect_true(result["dbovendry", "ph1to1h2o"] != 0)
  expect_true(all(eigen(result, only.values = TRUE, symmetric = TRUE)$values > 0))
})

test_that("build_kssl_fallback_matrix() includes real, non-identity ilr1/ilr2 correlations (the ILR-convention fix)", {
  result <- build_kssl_fallback_matrix(c("dbovendry", "ilr1", "ilr2"), genhz = "B")
  expect_true(result["dbovendry", "ilr1"] != 0)
  expect_true(result["dbovendry", "ilr2"] != 0)
  expect_true(result["ilr1", "ilr2"] != 0)
  expect_true(all(eigen(result, only.values = TRUE, symmetric = TRUE)$values > 0))
})

test_that("build_kssl_fallback_matrix() with genhz=NULL pools an unweighted average across genhz keys", {
  result <- build_kssl_fallback_matrix(c("dbovendry", "ph1to1h2o"), genhz = NULL)
  expect_equal(result["dbovendry", "ph1to1h2o"], 0.156698, tolerance = 1e-5)
})

test_that("build_kssl_fallback_matrix() returns NULL for a genhz absent from the property-matrix keys", {
  expect_null(build_kssl_fallback_matrix(c("dbovendry", "ph1to1h2o"), genhz = "R"))
  expect_null(build_kssl_fallback_matrix(c("dbovendry", "ph1to1h2o"), genhz = "nonexistent"))
})

test_that("build_kssl_fallback_matrix() returns NULL when fewer than 2 requested properties map to a KSSL column", {
  expect_null(build_kssl_fallback_matrix(c("dbovendry", "unmapped_a", "unmapped_b"), genhz = "A"))
  expect_null(build_kssl_fallback_matrix(character(0), genhz = "A"))
})

test_that("build_kssl_fallback_matrix() dimnames match the requested properties vector exactly, in order", {
  props <- c("rfv", "dbovendry", "cec7")
  result <- build_kssl_fallback_matrix(props, genhz = "C")
  expect_equal(rownames(result), props)
  expect_equal(colnames(result), props)
})

# ---------------------------------------------------------------------------
# 5 P2 chemistry properties (MULTI_PROPERTY_FUSION_PLAN.md task P2, KSSL-correlation design
# decision confirmed 2026-09-04: identity fallback via .kssl_property_name_map + this function).
# ---------------------------------------------------------------------------

test_that("build_kssl_fallback_matrix() treats the 5 P2 chemistry properties as identity (mapped in .kssl_property_name_map, absent from the real KSSL matrix)", {
  # Named in .kssl_property_name_map (unlike "unmapped_prop" above) but still not present in
  # kssl_property_matrices' own rownames - build_kssl_fallback_matrix()'s
  # `mapped_kssl_names %in% rownames(kssl_source)` filter must drop them the same way, not error.
  result <- build_kssl_fallback_matrix(c("dbovendry", "ph1to1h2o", "caco3", "ec", "ecec", "gypsum", "sar"), genhz = "A")
  for (p in c("caco3", "ec", "ecec", "gypsum", "sar")) {
    expect_equal(unname(result[p, p]), 1, info = p)
    expect_equal(unname(result["dbovendry", p]), 0, info = p)
  }
  expect_true(unname(result["dbovendry", "ph1to1h2o"]) != 0)  # real KSSL correlation unaffected
  expect_true(all(eigen(result, only.values = TRUE, symmetric = TRUE)$values > -1e-10))
})

test_that("build_kssl_fallback_matrix() stays positive-definite for every genhz with the full 14-property P2 vocabulary", {
  full_param_order <- c("db", "wr_3b", "wr_15b", "ilr1", "ilr2", "rfv", "ph", "cec", "soc",
                        "caco3", "ec", "ecec", "gypsum", "sar")
  for (g in c("O", "A", "E", "B", "C", "Cr")) {
    result <- build_kssl_fallback_matrix(full_param_order, genhz = g)
    expect_false(is.null(result), info = g)
    expect_equal(dim(result), c(14, 14), info = g)
    expect_true(all(eigen(result, only.values = TRUE, symmetric = TRUE)$values > -1e-10), info = g)
  }
})
