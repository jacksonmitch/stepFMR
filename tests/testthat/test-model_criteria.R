test_that("count_params_fmr counts parameters correctly", {
  expect_equal(
    count_params_fmr(ncol_het = 2, ncol_common = 1, G = 3, family = "gaussian"),
    3 * 2 + 1 + 2 + 3  # G*p_het + p_com + (G-1) + G sigmas
  )
  expect_equal(
    count_params_fmr(ncol_het = 2, ncol_common = 1, G = 3, family = "poisson"),
    3 * 2 + 1 + 2  # no sigma term
  )
})

test_that("compute_bic matches the standard formula", {
  expect_equal(compute_bic(loglik = -100, n = 50, k = 5), 200 + 5 * log(50))
})

test_that("select_best_G picks the lowest-BIC / highest-loglik fit", {
  fits <- list(
    list(bic = 120, loglik = -50),
    list(bic = 100, loglik = -55),
    list(bic = 150, loglik = -40)
  )
  expect_equal(select_best_G(fits, criterion = "bic")$bic, 100)
  expect_equal(select_best_G(fits, criterion = "loglik")$loglik, -40)
})