devtools::load_all()

match_groups <- function(tau, true_group) {
  n <- nrow(tau)
  g_hat <- ncol(tau)
  true_g <- max(true_group)

  indicator <- matrix(0, nrow = n, ncol = true_g)
  indicator[cbind(seq_len(n), true_group)] <- 1

  cost <- matrix(0, nrow = g_hat, ncol = true_g)
  for (k in seq_len(g_hat)) {
    for (t in seq_len(true_g)) {
      cost[k, t] <- sum((tau[, k] - indicator[, t])^2)
    }
  }

  if (g_hat <= true_g) {
    match_true <- as.integer(clue::solve_LSAP(cost))          # length g_hat
    match_fitted <- rep(NA_integer_, true_g)
    match_fitted[match_true] <- seq_len(g_hat)
  } else {
    match_fitted <- as.integer(clue::solve_LSAP(t(cost)))     # length true_g
  }

  match_fitted
}

run_replications <- function(params, n, control = fmr_control(), verbose = TRUE) {
  results <- vector("list", n_reps)

  for (i in seq_len(n_reps)) {
    if (verbose) cat(sprintf("Replication %d / %d...\n", i, n_reps))
    sim_data <- do.call(simulate_fmr, c(params, list(n = n, seed = i)))
    results[[i]] <- run_one_replication(
      rep_id = i, sim_data = sim_data, params = params, G_max = G_max, 
      control = control
    )
  }
 
  do.call(rbind, results)
}

run_one_replication <- function(rep_id, sim_data, params, G_max, control = fmr_control()) {
  t0 <- Sys.time()
  fit <- tryCatch(
    stepFMR(
      formula = sim_data$formula,
      data = sim_data$data,
      G_max = G_max,
      family = "gaussian",
      control = control
    ),
    error = function(e) e
  )
  runtime <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
 
  if (inherits(fit, "error")) {
    return(data.frame(
      rep = rep_id, runtime = runtime, error = TRUE,
      error_message = conditionMessage(fit),
      correct_selection = NA, correct_het = NA, correct_hom = NA,
      correct_outcome = NA, G_hat = NA,
      rmse_het = NA, rmse_hom = NA
    ))
  }
 
  selected <- fit$variable_selection$selected
  het_hat <- fit$effect_determination$heterogeneous
  hom_hat <- fit$effect_determination$homogeneous
 
  correct_selection <- isTRUE(setequal(selected, union(true_het, true_hom)))
  correct_het <- isTRUE(setequal(het_hat, true_het))
  correct_hom <- isTRUE(setequal(hom_hat, true_hom))
  correct_outcome <- correct_selection && correct_het && correct_hom
 
  true_G <- length(params$pi)
 
  rmse_het <- NA_real_
  rmse_hom <- NA_real_
  
  beta_g_hat <- fit$best_fit$parameter_values[["beta_g"]]
  beta_hat <- fit$best_fit$parameter_values[["beta"]]
  het_truth <- params$betas[, true_het, drop = FALSE]

  G_hat <- fit$best_fit$G
  tau <- fit$best_fit$parameter_values[["tau"]]

  match <- match_groups(tau, sim_data$true_group)

  # RMSE's
  errors <- list()

  for (v in true_het) {
    true_vals <- params$betas[, v]
    err <- if (!(v %in% selected)) {
      true_vals - 0
    } else if (v %in% het_hat) {
      fitted_vals <- rep(NA_real_, true_G)
      matched <- !is.na(match)
      fitted_vals[matched] <- beta_g_hat[match[matched], v]
      fitted_vals - true_vals
    } else {
      rep(beta_hat[[v]], true_G) - true_vals
    }
    errors[[v]] <- data.frame(variable = v, true_role = "het", sq_error = err^2)
  }

  for (v in true_hom) {
    true_val <- params$betas[1, v]
    err <- if (!(v %in% selected)) {
      true_val - 0
    } else if (v %in% hom_hat) {
      beta_hat[[v]] - true_val
    } else {
      beta_g_hat[, v] - true_val   # hom mistakenly split into het
    }
    errors[[v]] <- data.frame(variable = v, true_role = "hom", sq_error = err^2)
  }

  errors <- do.call(rbind, errors)

  rmse_het <- sqrt(mean(errors$sq_error[errors$true_role == "het"], na.rm = TRUE))
  rmse_hom <- sqrt(mean(errors$sq_error[errors$true_role == "hom"], na.rm = TRUE))

  data.frame(
    rep = rep_id, runtime = runtime, error = FALSE, error_message = NA,
    correct_selection = correct_selection, correct_het = correct_het,
    correct_hom = correct_hom, correct_outcome = correct_outcome,
    G_hat = G_hat, rmse_het = rmse_het, rmse_hom = rmse_hom
  )
}

summarize_results <- function(results) {
  ok <- results[!results$error, ]

  cat(sprintf("Replications: %d  |  Errors: %d\n\n", nrow(results), sum(results$error)))

  cat("--- Runtime ---\n")
  cat(sprintf("Min: %.2fs  |  Mean: %.2fs  |  Median: %.2fs  |  Max: %.2fs\n\n",
    min(ok$runtime), mean(ok$runtime), stats::median(ok$runtime), max(ok$runtime)
  ))

  cat("--- Model Selection ---\n")
  cat(sprintf("Completely Correct outcome:      %.1f%%\n", 100 * mean(ok$correct_outcome, na.rm = TRUE)))
  cat(sprintf("  Correct variable selection: %.1f%%\n", 100 * mean(ok$correct_selection, na.rm = TRUE)))
  cat(sprintf("  Correct heterogeneous set:  %.1f%%\n", 100 * mean(ok$correct_het, na.rm = TRUE)))
  cat(sprintf("  Correct homogeneous set:    %.1f%%\n", 100 * mean(ok$correct_hom, na.rm = TRUE)))
  cat(sprintf("G recovered correctly (G_hat == true G): %.1f%%\n\n",
    100 * mean(ok$G_hat == true_g, na.rm = TRUE)
  ))

  cat("--- Coefficient recovery (RMSE) ---\n")
  cat("Misclassifications and estimates compared against the true values\n")
  cat(sprintf("  RMSE heterogeneous coefficients: %.4f\n", mean(ok$rmse_het, na.rm = TRUE)))
  cat(sprintf("  RMSE homogeneous coefficients:   %.4f\n", mean(ok$rmse_hom, na.rm = TRUE)))
  
  invisible(ok)
}

true_het <- c("het1", "het2", "het3", "het4", "het5")
true_hom <- c("hom1", "hom2", "hom3")
true_null <- c("null1", "null2", "null3", "null4")

G_max <- 7
n_reps <- 100
params <- scenarios$four_group_twelve_variables_doubled
true_g <- length(params$pi)

# Remember to devtools::install() to reflect changes in parallel
NoSleepR::with_nosleep({
  # sim_results_500 <- run_replications(params, 500, verbose = FALSE)
  # cat("500 sequential\n")
  # summarize_results(sim_results_500)
  # sim_results_1000 <- run_replications( params, 1000, verbose = FALSE)
  # cat("1000 sequential\n")
  # summarize_results(sim_results_1000)
  sim_results_500_parallel <- run_replications(params, 
    500, fmr_control(parallel = TRUE), verbose = FALSE)
  cat("500 parallel\n")
  summarize_results(sim_results_500_parallel)
  # sim_results_1000_parallel <- run_replications(params,
  #    1000, fmr_control(parallel = TRUE), verbose = FALSE)
  # cat("1000 parallel\n")
  # summarize_results(sim_results_1000_parallel)
})

