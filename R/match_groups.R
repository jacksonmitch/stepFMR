#' Relabel a Fit's Mixture Components to Match Known True Groups
#'
#' Finds the best correspondence between a fitted model's mixture
#' components and a vector of known true group labels, then permutes the 
#' component-indexed parts of the fit so that they line up with the
#' true group ordering. 
#'
#' @param fit A `fit_fmr` object, or a top-level `stepFMR` result (in which
#'   case `fit$best_fit` is relabeled).
#' @param true_group An integer vector of true group labels, between 1:G,
#'  one entry per observation, as produced by `simulate_fmr()`.
#'
#' @return An object of the same class, with elements relabeled to true-group
#'   order.
#' @export
match_groups <- function(fit, true_group) {
  UseMethod("match_groups")
}

#' @export
match_groups.fit_fmr <- function(fit, true_group) {
  param_values <- fit$parameter_values
  tau <- param_values$tau

  checkmate::assert_matrix(tau, min.cols = 1)
  checkmate::assert_integer(true_group,
    len = nrow(tau), any.missing = FALSE, lower = 1
  )

  true_g <- max(true_group)

  g_hat <- ncol(tau)
  perm <- find_group_match(tau, true_group)
  matched <- !is.na(perm)

  g_names <- paste0("g", seq_len(true_g))

  new_tau <- matrix(NA_real_,
    nrow = nrow(tau), ncol = true_g,
    dimnames = list(NULL, g_names)
  )
  new_tau[, matched] <- tau[, perm[matched], drop = FALSE]
  param_values$tau <- new_tau

  new_pi_g <- stats::setNames(rep(NA_real_, true_g), g_names)
  new_pi_g[matched] <- param_values$pi_g[perm[matched]]
  param_values$pi_g <- new_pi_g

  new_beta_g <- matrix(NA_real_,
    nrow = true_g, ncol = ncol(param_values$beta_g),
    dimnames = list(g_names, colnames(param_values$beta_g))
  )
  new_beta_g[matched, ] <- param_values$beta_g[perm[matched], , drop = FALSE]
  param_values$beta_g <- new_beta_g

  if (!is.null(param_values$sigma_g)){
    new_sigma_g <- stats::setNames(rep(NA_real_, true_g), g_names)
    new_sigma_g[matched] <- param_values$sigma_g[perm[matched]]
    param_values$sigma_g <- new_sigma_g
  }
  
  fit$parameter_values <- param_values
  fit
}

#' @export
match_groups.stepFMR <- function(fit, true_group) {
  fit$best_fit <- match_groups(fit$best_fit, true_group)
  fit
}

#' @export
match_groups.default <- function(fit, true_group) {
  stop("match_groups() can't handle an object of class '",
    class(fit)[[1]], "'.", call. = FALSE
  )
}

find_group_match <- function(tau, true_group) {
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
    match_true <- as.integer(clue::solve_LSAP(cost))
    match_fitted <- rep(NA_integer_, true_g)
    match_fitted[match_true] <- seq_len(g_hat)
  } else {
    match_fitted <- as.integer(clue::solve_LSAP(t(cost)))
  }

  match_fitted
}