test_that("stepFMR validates inputs", {
  dat <- do.call(
    simulate_fmr,
    c(scenarios$two_group_effects, list(n = 200, seed = 1))
  )$data

  expect_error(stepFMR(y ~ x1, data = dat, G_max = 1))
  expect_error(stepFMR(y ~ x1, data = dat, family = "gamma"))
  expect_error(stepFMR(y ~ x1, data = dat, procedure = "both"))
  expect_error(stepFMR(y ~ x1, data = dat, control = list()))
  expect_error(stepFMR("not a formula", data = dat))
  expect_error(stepFMR(y ~ x1, data = as.matrix(dat)))
})

test_that("stepFMR runs end-to-end for each family", {
  gauss <- do.call(
    simulate_fmr,
    c(scenarios$two_group_effects, list(n = 300, seed = 1))
  )
  pois <- do.call(
    simulate_fmr,
    c(scenarios$two_group_effects_poisson, list(n = 300, seed = 2))
  )
  binom <- do.call(
    simulate_fmr,
    c(scenarios$two_group_effects_binomial, list(n = 300, seed = 3))
  )

  ctrl <- fmr_control(n_init = 5, n_kmeans_init = 2)

  fit_gauss <- stepFMR(gauss$formula, data = gauss$data, G_max = 3,
    family = "gaussian", control = ctrl)
  fit_pois <- stepFMR(pois$formula, data = pois$data, G_max = 3,
    family = "poisson", control = ctrl)
  fit_binom <- stepFMR(binom$formula, data = binom$data, G_max = 3,
    family = "binomial", control = ctrl)

  for (fit in list(fit_gauss, fit_pois, fit_binom)) {
    expect_s3_class(fit, "stepFMR")
    expect_s3_class(fit$best_fit, "fit_fmr")
    expect_s3_class(fit$variable_selection, "select_variables")
    expect_s3_class(fit$effect_determination, "determine_effects")
    expect_true(fit$best_fit$G %in% fit$G_values)
  }
})

test_that("stepFMR accepts cbind() binomial formulas as well as 0/1 responses", {
  binom <- do.call(
    simulate_fmr,
    c(scenarios$two_group_effects_binomial, list(n = 300, seed = 3))
  )

  dat_grouped <- binom$data
  dat_grouped$trials <- binom$data$size

  expect_no_error(
    stepFMR(cbind(y, trials - y) ~ x1 + x2 + x3, data = dat_grouped,
      family = "binomial", G_max = 3, control = fmr_control())
  )

  dat_bernoulli <- dat_grouped
  dat_bernoulli$y <- as.integer(dat_bernoulli$y > 0)

  expect_no_error(
    expect_warning(
      stepFMR(y ~ x1 + x2 + x3, data = dat_bernoulli,
        family = "binomial", G_max = 3, control = fmr_control())
    )
  )
})

test_that("stepFMR procedure argument controls which sub-results are populated", {
  sim <- do.call(
    simulate_fmr,
    c(scenarios$two_group_effects, list(n = 300, seed = 1))
  )

  fit_vars_only <- stepFMR(sim$formula, data = sim$data, G_max = 3,
    procedure = "variables")
  expect_s3_class(fit_vars_only$variable_selection, "select_variables")
  expect_null(fit_vars_only$effect_determination)

  fit_effects_only <- stepFMR(sim$formula, data = sim$data, G_max = 3,
    procedure = "effects")
  expect_null(fit_effects_only$variable_selection)
  expect_s3_class(fit_effects_only$effect_determination, "determine_effects")
})

test_that("stepFMR skips effect determination when no predictors are selected", {
  sim <- do.call(
    simulate_fmr,
    c(scenarios$two_group_effects, list(n = 300, seed = 1))
  )

  # alpha = 0 makes every forward-selection candidate ineligible
  fit <- stepFMR(sim$formula, data = sim$data, G_max = 3,
    control = fmr_control(alpha = 0))

  expect_length(fit$variable_selection$selected, 0)
  expect_null(fit$effect_determination)
  expect_s3_class(fit$best_fit, "fit_fmr")
})