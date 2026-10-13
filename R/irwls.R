#' Iteratively Reweighted Least Squares M-Step
#'
#' M-step for the Poisson and binomial families when the design has common
#' predictors. Starting from the current coefficients (or, if none exist, a
#' weighted fit to a transformed response), repeatedly forms the working
#' response and weights for the log/logit link, multiplies the weights by
#' the responsibilities \code{tau}, and re-solves with \code{wls_sqr()}
#' until the largest coefficient change is below \code{control$irwls_tol}
#' or \code{control$irwls_max_iter} is reached. \code{em_state} is updated
#' in place.
#'
#' @param dat A \code{FMRData} object.
#' @param em_state An \code{EmState} object.
#' @param control A \code{FMRControl} object.
#' @param family Either \code{"poisson"} or \code{"binomial"}.
#'
#' @return The updated \code{em_state}, with \code{beta_g}, \code{beta},
#'   \code{eta}, \code{irwls_iterations} and \code{irwls_converged} set.
#' @noRd
irwls_fmr <- function(dat, em_state, control,
                      family = c("poisson", "binomial")) {
  family <- match.arg(family)
  A <- dat$X_het
  B <- dat$X_com
  # stopifnot(!is.null(dat$B))
  y <- dat$y

  tau <- em_state$tau
  beta_g <- em_state$beta_g
  beta <- em_state$beta

  n <- nrow(A)
  G <- ncol(tau)

  maxit <- control$irwls_max_iter
  tol <- control$irwls_tol
  weight_floor <- control$weight_floor

  # I think trials should be saved somewhere but don't know where
  binomial_size <- dat$binomial_size
  # this will keep it being generated until we decide where it goes, but
  # calculations stay safe
  if (is.null(binomial_size)) {
    binomial_size <- rep(1, n)
  }

  binomial_size <- as.numeric(binomial_size)

  # Initial values for beta_g and beta
  if (is.null(beta_g) || is.null(beta)) {
    if (family == "poisson") {
      z0 <- log(pmax(y, 0) + 0.25)
    }

    if (family == "binomial") {
      p0 <- (y + 0.5) / (binomial_size + 1)
      p0 <- pmin(pmax(p0, 1e-8), 1 - 1e-8)
      z0 <- stats::qlogis(p0)
    }
    out <- wls_sqr(
      dat = dat,
      z = z0,
      weights = tau,
      control = control
    )
    beta_g <- out[["beta_g"]]
    beta <- out[["beta"]]
  }

  converged <- FALSE

  for (r in seq_len(maxit)) {
    beta_g_old <- beta_g
    beta_old <- beta

    eta <- linear_predictor_matrix(A, B, beta_g, beta)

    eta <- as.matrix(eta)
    storage.mode(eta) <- "double"

    if (family == "poisson") {
      eta_safe <- pmin(pmax(eta, -20), 20)

      mu <- exp(eta_safe)
      mu <- pmin(pmax(mu, 1e-10), 1e8)

      y_mat <- matrix(y, nrow = n, ncol = G)

      z <- eta_safe + (y_mat - mu) / mu
      w <- tau * mu
    }

    if (family == "binomial") {
      eta_safe <- pmin(pmax(eta, -30), 30)
      mu <- stats::plogis(eta_safe)
      mu <- pmin(pmax(mu, 1e-8), 1 - 1e-8)

      y_bar <- y / binomial_size
      y_bar_mat <- matrix(y_bar, nrow = n, ncol = G)

      z <- eta_safe + (y_bar_mat - mu) / pmax(mu * (1 - mu), 1e-10)
      w <- tau * binomial_size * mu * (1 - mu)
    }

    w <- pmax(w, weight_floor)


    out <- wls_sqr(
      dat = dat,
      z = z,
      weights = w,
      control = control
    )
    beta_g <- out[["beta_g"]]
    beta <- out[["beta"]]

    diff <- max(abs(c(beta_g - beta_g_old, beta - beta_old)))

    if (is.finite(diff) && diff < tol) {
      converged <- TRUE
      break
    }
  }

  em_state$beta_g <- beta_g
  em_state$beta <- beta
  em_state$eta <- eta
  em_state$irwls_iterations <- r
  em_state$irwls_converged <- converged
  em_state
}
