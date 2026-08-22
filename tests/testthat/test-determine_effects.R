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
    control = fmr_control()
  )

  expect_no_error(determine_effects(model = model, direction = "forward"))
  expect_no_error(determine_effects(model = model, direction = "backward"))
})

test_that("determine_effects correctly classifies effects (gaussian)", {
  sim <- do.call(
    simulate_fmr,
    c(scenarios$two_group_effects, list(n = 800, seed = 1))
  )

  model <- FMRModel$new(
    formula = sim$formula, data = sim$data,
    G_values = 2, family = "gaussian",
    control = fmr_control(n_init = 5, n_kmeans_init = 2)
  )

  fwd <- determine_effects(model = model, direction = "forward")
  bwd <- determine_effects(model = model, direction = "backward")

  # x1 differs by group in the true betas; x2 and x3 do not
  expect_true("x1" %in% fwd$heterogeneous)
  expect_true(all(c("x2", "x3") %in% fwd$homogeneous))
  expect_true("x1" %in% bwd$heterogeneous)
  expect_true(all(c("x2", "x3") %in% bwd$homogeneous))
})

test_that("determine_effects classifies effects correctly for poisson", {
  sim <- do.call(
    simulate_fmr,
    c(scenarios$two_group_effects_poisson, list(n = 800, seed = 2))
  )

  model <- FMRModel$new(
    formula = sim$formula, data = sim$data,
    G_values = 2, family = "poisson",
    control = fmr_control(n_init = 5, n_kmeans_init = 2)
  )

  fwd <- determine_effects(model = model, direction = "forward")
  expect_true("x1" %in% fwd$heterogeneous)
  expect_true(all(c("x2", "x3") %in% fwd$homogeneous))
})

test_that("determine_effects classifies effects correctly for binomial", {
  sim <- do.call(
    simulate_fmr,
    c(scenarios$two_group_effects_binomial, list(n = 800, seed = 3))
  )
  binom <- coerce_binomial_formula(sim$formula, sim$data)

  model <- FMRModel$new(
    formula = binom$formula, data = binom$data,
    G_values = 2, family = "binomial",
    control = fmr_control(n_init = 5, n_kmeans_init = 2)
  )

  fwd <- determine_effects(model = model, direction = "forward")
  expect_true("x1" %in% fwd$heterogeneous)
  expect_true(all(c("x2", "x3") %in% fwd$homogeneous))
})