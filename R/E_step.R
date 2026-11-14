#' E-Step
#'
#' Computes each component's mixing-weighted log-density for every observation
#' under the given family, uses log-sum-exp normalization to update the
#' responsibilities \code{tau}, and updates the mixing proportions
#' \code{pi_g} (floored at 1e-16 and renormalized). For \code{"zip"}, the
#' first \code{G - 1} components are Poisson and component \code{G} is a
#' point mass at zero. \code{em_state} is mutated in place.
#'
#' @param dat A \code{FMRData} object.
#' @param em_state An \code{EmState} object holding \code{pi_g},
#'   \code{sigma_g} (Gaussian only) and the linear predictor \code{eta}.
#' @param family One of \code{"gaussian"}, \code{"poisson"},
#'   \code{"binomial"}, or \code{"zip"}.
#'
#' @return Called for its side effect of updating \code{em_state}
#'   (\code{loglik}, \code{tau}, \code{pi_g}, \code{eta}).
#' @noRd
e_step <- function(dat, em_state, family) {
  y <- dat$y
  n <- dat$n

  pi_g <- em_state$pi_g
  sigma_g <- em_state$sigma_g
  G <- em_state$G
  eta <- em_state$eta

  log_pi <- log(pi_g)
  log_weights <- matrix(NA_real_, nrow = n, ncol = G)

  # Family specific E-step quantities

  if (family == "gaussian") {
    stopifnot(length(sigma_g) == G)

    for (g in seq_len(G)) {
      log_weights[, g] <- log_pi[g] + stats::dnorm(
        x = y,
        mean = eta[, g],
        sd = pmax(sigma_g[g], 1e-16),
        log = TRUE
      )
    }
  } else if (family == "poisson") {
    mu <- exp(eta)

    for (g in seq_len(G)) {
      log_weights[, g] <- log_pi[g] + stats::dpois(
        x = y,
        lambda = mu[, g],
        log = TRUE
      )
    }
  } else if (family == "binomial") {
    binomial_size <- dat$binomial_size

    mu <- stats::plogis(eta)
    mu <- pmin(pmax(mu, 1e-8), 1 - 1e-8)

    for (g in seq_len(G)) {
      log_weights[, g] <- log_pi[g] + stats::dbinom(
        x = y,
        size = binomial_size,
        prob = mu[, g],
        log = TRUE
      )
    }
  } else if (family == "zip") {
    mu <- exp(eta)

    for (g in seq_len(G-1)) {
      log_weights[, g] <- log_pi[g] + stats::dpois(
        x = y,
        lambda = mu[, g],
        log = TRUE
      )
    }
    log_weights[, G] <- log_pi[G] + ifelse(y == 0, 0, -Inf)
  }

  # Shift each row before exponentiating so small densities do not underflow.
  log_weight_max <- log_weights[, 1L]
  if (G > 1L) {
    for (g in 2:G) {
      log_weight_max <- pmax(log_weight_max, log_weights[, g])
    }
  }
  if (any(!is.finite(log_weight_max))) {
    stop(sprintf(
      "%s mixture has no finite log-density for one or more observations.",family
    ), call. = FALSE)
  }

  weights <- exp(log_weights - log_weight_max)
  sum_weights <- rowSums(weights)
  loglik <- sum(log_weight_max + log(sum_weights))
  tau <- weights / sum_weights

  pi_new <- colMeans(tau)
  pi_new <- pmax(pi_new, 1e-16)
  pi_new <- pi_new / sum(pi_new)

  em_state$loglik <- loglik
  em_state$tau <- tau
  em_state$pi_g <- pi_new
  em_state$eta <- eta
}
