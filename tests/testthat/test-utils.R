test_that("coerce_binomial_formula handles cbind(success, failure) formulas", {
  dat <- data.frame(success = c(3, 5), failure = c(7, 5), x = c(1, 2))
  out <- coerce_binomial_formula(cbind(success, failure) ~ x, dat)

  expect_equal(out$data$.binom_size, c(10, 10))
  expect_identical(deparse(out$formula[[2]]), "success")
})

test_that("coerce_binomial_formula handles plain 0/1 responses", {
  dat <- data.frame(y = c(0, 1, 1, 0), x = c(1, 2, 3, 4))
  out <- coerce_binomial_formula(y ~ x, dat)

  expect_equal(out$data$.binom_size, rep(1, 4))
})

test_that("coerce_binomial_formula errors on invalid inputs", {
  expect_error(
    coerce_binomial_formula(y ~ x, data.frame(y = c(0, 1, 2), x = c(1, 2, 3)))
  )
  expect_error(
    coerce_binomial_formula(
      y ~ x, data.frame(y = c(0, 1), x = c(1, 2), .binom_size = c(1, 1))
    )
  )
  expect_error(
    coerce_binomial_formula(
      cbind(success, failure) ~ 1,
      data.frame(success = c(-1, 2), failure = c(3, 4))
    )
  )
})