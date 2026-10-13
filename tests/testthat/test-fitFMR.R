ctrl <- fmr_control(n_init = 5, n_kmeans_init = 2)

gauss <- do.call(
  simulate_fmr,
  c(scenarios$two_group_effects, list(n = 300, seed = 1))
)

test_that("fitFMR validates inputs", {
  expect_error(fitFMR(gauss$formula, gauss$data, G = 1), "G")
  expect_error(fitFMR(gauss$formula, gauss$data, G = 2.5))
  expect_error(fitFMR(gauss$formula, gauss$data, family = "gamma"), "family")
  expect_error(fitFMR(gauss$formula, gauss$data, family = "zip"), "family")
  expect_error(fitFMR(gauss$formula, gauss$data, common = "zz"), "zz")
  expect_error(
    fitFMR(gauss$formula, gauss$data, control = list()),
    "control"
  )
  expect_error(fitFMR("not a formula", gauss$data))
  expect_error(fitFMR(gauss$formula, as.matrix(gauss$data)))
})

test_that("fitFMR returns the documented object for a vector of G", {
  fit <- fitFMR(gauss$formula, gauss$data, G = 2:3, control = ctrl)

  expect_s3_class(fit, "fitFMR")
  expect_named(fit, c(
    "call", "formula", "family", "G_values", "common", "heterogeneous",
    "fits", "table", "selected_G", "best_fit", "control"
  ))
  expect_equal(fit$G_values, 2:3)
  expect_named(fit$fits, c("G2", "G3"))
  expect_true(all(vapply(fit$fits, inherits, logical(1), what = "fit_fmr")))
  expect_equal(
    names(fit$table),
    c("G", "loglik", "num_parameters", "bic", "iterations", "converged")
  )
  expect_equal(fit$table$G, 2:3)
  expect_equal(fit$selected_G, fit$table$G[which.min(fit$table$bic)])
  expect_identical(
    fit$best_fit,
    fit$fits[[paste0("G", fit$selected_G)]]
  )
  expect_equal(fit$heterogeneous, c("x1", "x2", "x3"))
  expect_equal(fit$common, character(0))
})

test_that("fitFMR with a scalar G fits exactly one model", {
  fit <- fitFMR(gauss$formula, gauss$data, G = 2, control = ctrl)

  expect_length(fit$fits, 1L)
  expect_equal(nrow(fit$table), 1L)
  expect_equal(fit$selected_G, 2L)
})

test_that("fitFMR sorts and de-duplicates G", {
  fit <- fitFMR(gauss$formula, gauss$data, G = c(3, 2, 3), control = ctrl)

  expect_equal(fit$G_values, 2:3)
  expect_named(fit$fits, c("G2", "G3"))
})

test_that("common predictors get one shared coefficient", {
  fit <- fitFMR(
    gauss$formula, gauss$data,
    G = 2, common = "x2", control = ctrl
  )
  pv <- fit$best_fit$parameter_values

  expect_equal(fit$common, "x2")
  expect_equal(fit$heterogeneous, c("x1", "x3"))
  expect_equal(names(pv$beta), "x2")
  expect_equal(colnames(pv$beta_g), c("(Intercept)", "x1", "x3"))
})

test_that("fitFMR works when every predictor is common", {
  fit <- fitFMR(
    gauss$formula, gauss$data,
    G = 2, common = c("x1", "x2", "x3"), control = ctrl
  )
  pv <- fit$best_fit$parameter_values

  expect_equal(fit$heterogeneous, character(0))
  expect_equal(colnames(pv$beta_g), "(Intercept)")
  expect_equal(names(pv$beta), c("x1", "x2", "x3"))
})

test_that("fitFMR runs for poisson and binomial responses", {
  pois <- do.call(
    simulate_fmr,
    c(scenarios$two_group_effects_poisson, list(n = 300, seed = 2))
  )
  binom <- do.call(
    simulate_fmr,
    c(scenarios$two_group_effects_binomial, list(n = 300, seed = 3))
  )

  fit_pois <- fitFMR(pois$formula, pois$data,
    G = 2:3, family = "poisson", common = "x2", control = ctrl
  )
  fit_binom <- fitFMR(binom$formula, binom$data,
    G = 2:3, family = "binomial", common = "x2", control = ctrl
  )

  for (fit in list(fit_pois, fit_binom)) {
    expect_s3_class(fit, "fitFMR")
    expect_true(fit$selected_G %in% 2:3)
    expect_s3_class(fit$best_fit, "fit_fmr")
  }
  expect_equal(fit_pois$family, "poisson")
  expect_equal(fit_binom$family, "binomial")
})

test_that("fitFMR warns when EM does not converge", {
  # a tiny tol with only 1 iteration cannot meet the relative-change
  # stopping rule
  no_conv <- fmr_control(
    n_init = 3, n_kmeans_init = 1, tol = 1e-16, max_iter = 1
  )

  expect_warning(
    fitFMR(gauss$formula, gauss$data, G = 2, control = no_conv),
    "did not converge"
  )
})

test_that("print.fitFMR shows the BIC table and the selected fit", {
  fit <- fitFMR(
    gauss$formula, gauss$data,
    G = 2:3, common = "x2", control = ctrl
  )
  output <- capture.output(res <- print(fit))

  expect_identical(res, fit)
  expect_match(output[[1]], "^fitFMR \\| family: gaussian \\| G: 2,3")
  expect_true(any(grepl("Heterogeneous: x1, x3", output)))
  expect_true(any(grepl("Homogeneous:   x2", output)))
  expect_true(any(grepl("bic", output)))
  expect_equal(sum(grepl("\\*", output)), 1L)
  expect_true(any(grepl("Best fit", output)))
})

test_that("summary.fitFMR collects the results and prints them", {
  fit <- fitFMR(
    gauss$formula, gauss$data,
    G = 2:3, common = "x2", control = ctrl
  )
  s <- summary(fit)

  expect_s3_class(s, "summary.fitFMR")
  expect_named(s, c(
    "call", "formula", "family", "G_values", "common", "heterogeneous",
    "table", "selected_G", "best_fit", "initialization"
  ))
  expect_equal(s$selected_G, fit$selected_G)

  output <- capture.output(res <- print(s))
  expect_identical(res, s)
  expect_match(output[[1]], "^fitFMR summary \\| family: gaussian")
  expect_true(any(grepl("converged", output)))
  expect_true(any(grepl("Initialization \\(G = ", output)))
  expect_true(any(grepl("Best fit", output)))
})

test_that("print.summary.fitFMR notes non-convergence", {
  no_conv <- fmr_control(
    n_init = 3, n_kmeans_init = 1, tol = 1e-16, max_iter = 1
  )
  fit <- suppressWarnings(
    fitFMR(gauss$formula, gauss$data, G = 2, control = no_conv)
  )
  output <- capture.output(print(summary(fit)))

  expect_true(any(grepl("did not converge for G = 2", output)))
})

test_that("match_groups relabels the selected fit of a fitFMR object", {
  fit <- fitFMR(gauss$formula, gauss$data, G = 2, control = ctrl)
  aligned <- match_groups(fit, gauss$true_group)

  expect_s3_class(aligned, "fitFMR")
  expect_identical(aligned$fits, fit$fits)

  # scenario mixing proportions are c(0.4, 0.6): after matching, g1 is smaller
  pi_hat <- aligned$best_fit$parameter_values$pi_g
  expect_equal(names(pi_hat), c("g1", "g2"))
  expect_lt(pi_hat[["g1"]], pi_hat[["g2"]])
})
