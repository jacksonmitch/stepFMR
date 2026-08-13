test_data <- do.call(
  simulate_fmr,
  c(
    scenarios$two_group_effects,
    list(n = 400, seed = 123)
  )
)$data

test_that("stepFMR runs without error with control$parallel = TRUE", {
  expect_no_error(
    stepFMR(formula(test_data), test_data,
      control = build_control(parallel = TRUE)
    )
  )
})

test_that("stepFMR runs without error with user specified plan", {
  future::plan(future::multisession)
  expect_no_error(
    stepFMR(formula(test_data), test_data)
  )
  future::plan(future::sequential)
})
