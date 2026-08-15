test_that("gaussian fits accurately", {
  sim <- do.call(
    simulate_fmr,
    c(scenarios$two_group_effects, list(n = 1000, seed = 1))
  )

  model <- FMRModel$new(
    formula = sim$formula,
    data = sim$data,
    G_values = 2,
    family = "gaussian",
    control = build_control(n_init = 5, n_kmeans_init = 2, max_iter = 200)
  )

  prepared <- prepare_data(model, common = c("x2", "x3"))
  fit <- fit_across_G(model, prepared)[[1]]

  expect_true(fit$converged)

  fit <- match_groups(fit, sim$true_group)
  em_values <- fit$parameter_values
  true_betas <- scenarios$two_group_effects$betas

  het_truth <- unname(true_betas[, c("(Intercept)", "x1")])
  print(em_values$beta_g)
  fitted_het <- unname(em_values$beta_g[, c("(Intercept)", "x1")])

  fitted_common <- unname(em_values$beta)
  common_truth <- unname(true_betas["g1", c("x2", "x3")])

  expect_equal(round(fitted_het, 1), round(het_truth, 1),
    tolerance = 0.005
  )
  expect_equal(round(fitted_common, 1), round(common_truth, 1),
    tolerance = 0.005
  )
})

test_that("poisson fits accurately", {
  sim <- do.call(
    simulate_fmr,
    c(scenarios$two_group_effects_poisson, list(n = 1500, seed = 2))
  )

  model <- FMRModel$new(
    formula = sim$formula,
    data = sim$data,
    G_values = 2,
    family = "poisson",
    control = build_control(n_init = 5, n_kmeans_init = 2, max_iter = 200)
  )

  prepared <- prepare_data(model, common = c("x2", "x3"))
  fit <- fit_across_G(model, prepared)[[1]]

  expect_true(fit$converged)

  fit <- match_groups(fit, sim$true_group)
  em_values <- fit$parameter_values
  true_betas <- scenarios$two_group_effects_poisson$betas

  het_truth <- unname(true_betas[, c("(Intercept)", "x1")])
  fitted_het <- unname(em_values$beta_g[, c("(Intercept)", "x1")])

  fitted_common <- unname(em_values$beta)
  common_truth <- unname(true_betas["g1", c("x2", "x3")])

  expect_equal(fitted_het, het_truth, tolerance = 0.1)
  expect_equal(fitted_common, common_truth, tolerance = 0.1)
})

test_that("binomial fits accurately", {
  sim <- do.call(
    simulate_fmr,
    c(scenarios$two_group_effects_binomial, list(n = 1500, seed = 3))
  )
  fixed <- coerce_binomial_formula(sim$formula, sim$data)

  model <- FMRModel$new(
    formula = fixed$formula,
    data = fixed$data,
    G_values = 2,
    family = "binomial",
    control = build_control(n_init = 5, n_kmeans_init = 2, max_iter = 200)
  )

  prepared <- prepare_data(model, common = c("x2", "x3"))
  fit <- fit_across_G(model, prepared)[[1]]

  expect_true(fit$converged)

  fit <- match_groups(fit, sim$true_group)
  em_values <- fit$parameter_values
  true_betas <- scenarios$two_group_effects_binomial$betas

  het_truth <- unname(true_betas[, c("(Intercept)", "x1")])
  fitted_het <- unname(em_values$beta_g[, c("(Intercept)", "x1")])

  fitted_common <- unname(em_values$beta)
  common_truth <- unname(true_betas["g1", c("x2", "x3")])

  expect_equal(fitted_het, het_truth, tolerance = 0.02)
  expect_equal(fitted_common, common_truth, tolerance = 0.02)
})
