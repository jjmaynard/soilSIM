test_that("build_cache_key() produces a stable, filename-safe key for the same AOI/id/depth/kind", {
  aoi <- terra::vect("POLYGON((0 0, 1 0, 1 1, 0 1, 0 0))")
  key1 <- build_cache_key(aoi, "ph", 0, 20, "ssurgo")
  key2 <- build_cache_key(aoi, "ph", 0, 20, "ssurgo")
  expect_equal(key1, key2)
  expect_match(key1, "^[A-Za-z0-9_-]+$")
})

test_that("build_cache_key() differs when any input differs", {
  aoi1 <- terra::vect("POLYGON((0 0, 1 0, 1 1, 0 1, 0 0))")
  aoi2 <- terra::vect("POLYGON((2 2, 3 2, 3 3, 2 3, 2 2))")

  base_key <- build_cache_key(aoi1, "ph", 0, 20, "ssurgo")
  expect_false(base_key == build_cache_key(aoi2, "ph", 0, 20, "ssurgo"))
  expect_false(base_key == build_cache_key(aoi1, "clay", 0, 20, "ssurgo"))
  expect_false(base_key == build_cache_key(aoi1, "ph", 20, 50, "ssurgo"))
  expect_false(base_key == build_cache_key(aoi1, "ph", 0, 20, "solus"))
})

test_that("cache_set()/cache_get() round-trip a value and respect TTL", {
  aoi <- terra::vect("POLYGON((0 0, 1 0, 1 1, 0 1, 0 0))")
  key <- build_cache_key(aoi, paste0("test-", as.integer(Sys.time())), 0, 20, "unittest")
  on.exit(unlink(file.path(tools::R_user_dir("soilSIM", "cache"), paste0(key, ".rds"))), add = TRUE)

  expect_null(cache_get(key))

  value <- data.frame(x = 1:3)
  expect_true(cache_set(key, "unittest", value))
  expect_equal(cache_get(key), value)

  # A TTL of 0 seconds means anything just written is already "too old".
  expect_null(cache_get(key, ttl_seconds = 0))
})

test_that("cache_get() returns NULL (not an error) for a nonexistent key", {
  expect_null(cache_get("this-key-should-never-exist-12345"))
})
