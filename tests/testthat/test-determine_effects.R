test_that("determine_effects", {
  test_data <- do.call(
    simulate_fmr,
    c(
      scenarios$two_group_effects,
      list(n = 400, seed = 123)
    )
  )$data

  model <- FMRModel$new(
    formula = formula(test_data),
    data = test_data,
    G_values = 2,
    family = "gaussian",
    control = build_control()
  )

  expect_no_error(determine_effects(model = model, direction = "forward"))
  expect_no_error(determine_effects(model = model, direction = "backward"))
})
