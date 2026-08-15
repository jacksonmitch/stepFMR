test_that("determine_effects", {
  sim <- do.call(
    simulate_fmr,
    c(
      scenarios$two_group_effects,
      list(n = 400, seed = 123)
    )
  )

  model <- FMRModel$new(
    formula = sim$formula,
    data = sim$data,
    G_values = 2,
    family = "gaussian",
    control = build_control()
  )

  expect_no_error(determine_effects(model = model, direction = "forward"))
  expect_no_error(determine_effects(model = model, direction = "backward"))
})

test_that("determine_effects", {
  sim <- do.call(
    simulate_fmr,
    c(
      scenarios$two_group_effects,
      list(n = 400, seed = 123)
    )
  )

  model <- FMRModel$new(
    formula = sim$formula,
    data = sim$data,
    G_values = 2,
    family = "gaussian",
    control = build_control()
  )

  expect_no_error(determine_effects(model = model, direction = "forward"))
  expect_no_error(determine_effects(model = model, direction = "backward"))
})
