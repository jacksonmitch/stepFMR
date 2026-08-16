test_that("select_variables validates direction", {
  dat <- do.call(
    simulate_fmr,
    c(scenarios$two_group_effects, list(n = 100, seed = 1))
  )$data

  model <- FMRModel$new(
    formula = y ~ x1, data = dat,
    G_values = 2, family = "gaussian", control = build_control()
  )

  expect_error(select_variables(model = model, direction = "sideways"))
})

test_that("select_variables drops a null predictor (gaussian)", {
  sim <- do.call(
    simulate_fmr,
    c(scenarios$two_group_four_variables, list(n = 200, seed = 1))
  )

  model <- FMRModel$new(
    formula = sim$formula, data = sim$data,
    G_values = 2, family = "gaussian", control = build_control()
  )

  fwd <- select_variables(model = model, direction = "forward")
  expect_s3_class(fwd, "select_variables")
  expect_true(all(c("x1", "x2", "x3") %in% fwd$selected))
  expect_false("x4" %in% fwd$selected)

  bwd <- select_variables(model = model, direction = "backward")
  expect_true(all(c("x1", "x2", "x3") %in% bwd$selected))
  expect_false("x4" %in% bwd$selected)
})

test_that("select_variables drops a null predictor (poisson)", {
  sim <- do.call(
    simulate_fmr,
    c(scenarios$two_group_four_variables, list(n = 100, seed = 1,
    family = "poisson"))
  )

  model <- FMRModel$new(
    formula = sim$formula, data = sim$data,
    G_values = 2, family = "poisson", control = build_control()
  )

  fwd <- select_variables(model = model, direction = "forward")
  expect_true(all(c("x1", "x2", "x3") %in% fwd$selected))
  expect_false("x4" %in% fwd$selected)
})

test_that("select_variables runs without error for binomial", {
  sim <- do.call(
    simulate_fmr,
    c(scenarios$two_group_effects_binomial, list(n = 400, seed = 1))
  )
  binom <- coerce_binomial_formula(sim$formula, sim$data)

  model <- FMRModel$new(
    formula = binom$formula, data = binom$data,
    G_values = 2, family = "binomial",
    control = build_control()
  )

  expect_no_error(select_variables(model = model, direction = "forward"))
  expect_no_error(select_variables(model = model, direction = "backward"))
})