test_that("is_unsuitable() vectorized toupper()/trimws() matches the original per-row computation", {
  # Regression test for the PERFORMANCE_IMPROVEMENT_PLAN.md Tier 2 is_unsuitable() fix:
  # toupper()/trimws() are now computed once (vectorized) instead of per-row inside the
  # classification loop - this must be a pure vectorization with identical results.
  data <- data.frame(
    hzname = c(" a ", "Bw", "R", " r ", "  ", NA, "2Cr1", "Om", "Bt"),
    desgnmaster = c("A", "B", "R", "C", NA, "O", "C", NA, "B"),
    stringsAsFactors = FALSE
  )
  result <- is_unsuitable(data)
  expect_length(result, nrow(data))
  expect_true(result[3])   # "R" hzname
  expect_true(result[4])   # " r " trims/uppercases to "R", matches R$ pattern
  expect_true(result[6])   # NA hzname but desgnmaster "O"
  expect_false(result[2])  # "Bw" - ordinary suitable horizon
})

test_that("is_unsuitable() works when desgnmaster is absent (no crash, master treated as NA)", {
  data <- data.frame(hzname = c("A", "R", "Bw"), stringsAsFactors = FALSE)
  result <- is_unsuitable(data)
  expect_length(result, 3)
  expect_true(result[2])
  expect_false(result[1])
})

test_that("is_unsuitable() returns all FALSE when hzname_col is missing", {
  data <- data.frame(other_col = 1:3)
  expect_equal(is_unsuitable(data), rep(FALSE, 3))
})
