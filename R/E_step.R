# E-step for Gaussian mixture regression

e_step <- function(dat, em_state, family) {
  y <- dat$y
  n <- dat$n

  pi_g <- em_state$pi_g
  sigma_g <- em_state$sigma_g
  G <- em_state$G
  eta <- em_state$eta

  weights <- matrix(NA_real_, nrow = n, ncol = G)

  # Family specific E-step quantities

  if (family == "gaussian") {
    stopifnot(length(sigma_g) == G)

    for (g in seq_len(G)) {
      weights[, g] <- pi_g[g] * stats::dnorm(
        x = y,
        mean = eta[, g],
        sd = pmax(sigma_g[g], 1e-16),
        log = FALSE
      )
    }
  } else if (family == "poisson") {
    mu <- exp(eta)

    for (g in seq_len(G)) {
      weights[, g] <- pi_g[g] * stats::dpois(
        x = y,
        lambda = mu[, g],
        log = FALSE
      )
    }
  } else if (family == "binomial") {
    binomial_size <- dat$binomial_size

    mu <- stats::plogis(eta)
    mu <- pmin(pmax(mu, 1e-8), 1 - 1e-8)

    for (g in seq_len(G)) {
      weights[, g] <- pi_g[g] * stats::dbinom(
        x = y,
        size = binomial_size,
        prob = mu[, g],
        log = FALSE
      )
    }
  } else if (family == "zip") {
    mu <- exp(eta)

    for (g in seq_len(G-1)) {
      weights[, g] <- pi_g[g] * stats::dpois(
        x = y,
        lambda = mu[, g],
        log = FALSE
      )
    }
    weights[, G] <- pi_g[G] * ifelse(y == 0, 1, 0)
  }
  stopifnot(all(weights >= 0))

  sum_weights <- rowSums(weights)
  loglik <- sum(log(sum_weights))
  tau <- weights / sum_weights

  pi_new <- colMeans(tau)
  pi_new <- pmax(pi_new, 1e-16)
  pi_new <- pi_new / sum(pi_new)

  em_state$loglik <- loglik
  em_state$tau <- tau
  em_state$pi_g <- pi_new
  em_state$eta <- eta
}
