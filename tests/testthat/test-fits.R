test_that("gaussian fits accurately", {
  dat <- do.call(
    simulate_fmr,
    c(scenarios$two_group_effects, list(n = 1000, seed = 1))
  )$data

  model <- FMRModel$new(
    formula = y ~ x1 + x2 + x3,
    data = dat,
    G_values = 2,
    family = "gaussian",
    control = build_control(n_init = 5, n_kmeans_init = 2, max_iter = 200)
  )

  prepared <- prepare_data(model, common = c("x2", "x3"))
  fit <- fit_across_G(model, prepared)[[1]]

  expect_true(fit$converged)

  em_values <- fit$parameter_values
  true_betas <- scenarios$two_group_effects$betas

  het_truth <- unname(true_betas[, c("Intercept", "x1")])
  perm <- match_group_order(em_values$beta_g, het_truth)

  fitted_het <- unname(em_values$beta_g[perm, , drop = FALSE])
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
  dat <- do.call(
    simulate_fmr,
    c(scenarios$two_group_effects_poisson, list(n = 1500, seed = 2))
  )$data

  model <- FMRModel$new(
    formula = y ~ x1 + x2 + x3,
    data = dat,
    G_values = 2,
    family = "poisson",
    control = build_control(n_init = 5, n_kmeans_init = 2, max_iter = 200)
  )

  prepared <- prepare_data(model, common = c("x2", "x3"))
  fit <- fit_across_G(model, prepared)[[1]]

  expect_true(fit$converged)

  em_values <- fit$parameter_values
  true_betas <- scenarios$two_group_effects_poisson$betas

  het_truth <- unname(true_betas[, c("Intercept", "x1")])
  perm <- match_group_order(em_values$beta_g, het_truth)

  fitted_het <- unname(em_values$beta_g[perm, , drop = FALSE])
  fitted_common <- unname(em_values$beta)
  common_truth <- unname(true_betas["g1", c("x2", "x3")])

  expect_equal(fitted_het, het_truth, tolerance = 0.1)
  expect_equal(fitted_common, common_truth, tolerance = 0.1)
})

test_that("binomial fits accurately", {
  dat <- do.call(
    simulate_fmr,
    c(scenarios$two_group_effects_binomial, list(n = 1500, seed = 3))
  )
  dat$data$.binom_size <- dat$size

  model <- FMRModel$new(
    formula = y ~ x1 + x2 + x3,
    data = dat$data,
    G_values = 2,
    family = "binomial",
    control = build_control(n_init = 5, n_kmeans_init = 2, max_iter = 200)
  )

  prepared <- prepare_data(model, common = c("x2", "x3"))
  fit <- fit_across_G(model, prepared)[[1]]

  expect_true(fit$converged)

  em_values <- fit$parameter_values
  true_betas <- scenarios$two_group_effects_binomial$betas

  het_truth <- unname(true_betas[, c("Intercept", "x1")])
  perm <- match_group_order(em_values$beta_g, het_truth)

  fitted_het <- unname(em_values$beta_g[perm, , drop = FALSE])
  fitted_common <- unname(em_values$beta)
  common_truth <- unname(true_betas["g1", c("x2", "x3")])

  expect_equal(fitted_het, het_truth, tolerance = 0.02)
  expect_equal(fitted_common, common_truth, tolerance = 0.02)
})
