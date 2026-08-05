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
