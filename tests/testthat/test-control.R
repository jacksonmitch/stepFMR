test_that("build_control returns valid defaults", {
  ctrl <- build_control()
  expect_s3_class(ctrl, "FMRControl")
  expect_identical(ctrl$alpha, 0.05)
})

test_that("build_control rejects invalid alpha", {
  expect_error(build_control(alpha = 1.5))
  expect_error(build_control(alpha = -0.1))
})

test_that("build_control rejects n_best_init > n_init", {
  expect_error(build_control(n_init = 2, n_best_init = 5))
})

test_that("build_control rejects n_kmeans_init > n_init", {
  expect_error(build_control(n_init = 2, n_kmeans_init = 5))
})

test_that("build_control rejects invalid direction", {
  expect_error(build_control(direction = "sideways"))
})

test_that("build_control rejects negative sigma_floor", {
  expect_error(build_control(sigma_floor = -1))
})

test_that("build_control rejects a non-flag parallel argument", {
  expect_error(build_control(parallel = "yes"))
})