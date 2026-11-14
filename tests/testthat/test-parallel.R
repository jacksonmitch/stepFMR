sim <- do.call(
  simulate_fmr,
  c(
    scenarios$two_group_effects,
    list(n = 400, seed = 123)
  )
)

test_that("stepFMR runs without error with control$parallel = TRUE", {
  expect_no_error(
    stepFMR(sim$formula, sim$data,
      control = fmr_control(parallel = TRUE, verbose = FALSE)
    )
  )
})

test_that("stepFMR runs without error with user specified plan", {
  old_plan <- future::plan(future::multisession)
  on.exit(future::plan(old_plan), add = TRUE)

  expect_no_error(
    stepFMR(sim$formula, sim$data, control = fmr_control(verbose = FALSE))
  )
})
