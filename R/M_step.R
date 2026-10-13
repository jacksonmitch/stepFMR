# Wrapper function for M-Steps across families and method
#' M-Step Dispatcher
#'
#' Chooses the M-step for the family. If the design has common
#' (homogeneous) predictors (\code{dat$p_com > 0}), the structured-QR
#' version (\code{m_step_sqr_*}) is used. Otherwise the plain
#' per-component fit (\code{lm.wfit}/\code{glm.fit}) is used. The
#' \code{"zip"} family reuses the Poisson M-steps. \code{em_state} is
#' mutated in place.
#'
#' @param dat A \code{FMRData} object.
#' @param em_state An \code{EmState} object.
#' @param family One of \code{"gaussian"}, \code{"poisson"},
#'   \code{"binomial"}, or \code{"zip"}.
#' @param control A \code{FMRControl} object.
#'
#' @return Called for its side effect of updating \code{em_state}
#'   (\code{beta_g}, \code{beta}, \code{eta}, and \code{sigma_g} for
#'   Gaussian).
#' @noRd
m_step <- function(dat, em_state, family, control) {
  if (dat$p_com > 0) {
    em_state <- switch(family,
      gaussian = m_step_sqr_gaussian(dat, em_state, control),
      poisson  = m_step_sqr_poisson(dat, em_state, control),
      zip      = m_step_sqr_poisson(dat, em_state, control),
      binomial = m_step_sqr_binomial(dat, em_state, control),
      stop(sprintf("No sQR M-step for family '%s'", family))
    )
  } else {
    em_state <- switch(family,
      gaussian = m_step_gaussian(dat, em_state, control),
      poisson  = m_step_poisson(dat, em_state, control),
      zip      = m_step_poisson(dat, em_state, control),
      binomial = m_step_binomial(dat, em_state, control),
      stop(sprintf("No M-step for family '%s'", family))
    )
  }
}

# Gaussian M-step, no common predictors: one weighted least squares fit per
# component (lm.wfit, weights = tau), then component SDs from the weighted
# residuals, floored at control$sigma_floor.
m_step_gaussian <- function(dat, em_state, control) {
  tau <- em_state$tau
  G <- em_state$G
  X <- dat$X_het
  y <- dat$y
  p <- dat$p_het

  beta_g <- matrix(NA_real_, nrow = G, ncol = p)

  for (g in seq_len(G)) {
    wg <- pmax(tau[, g], 0)

    fit_g <- stats::lm.wfit(
      x = X,
      y = y,
      w = wg
    )

    coef_g <- fit_g$coefficients
    coef_g[is.na(coef_g)] <- 0

    beta_g[g, ] <- coef_g
  }

  # TODO: look at this line and think of alternatives
  mu <- linear_predictor_matrix(
    A = dat$X_het,
    B = dat$X_com,
    beta_g = beta_g,
    beta = NULL
  )
  res2 <- (dat$y - mu)^2

  # Would it be faster to do these with matrix form instead of col sums?
  # TODO: test with scripts (maybe do this when in C++ land I imagine)

  den <- pmax(colSums(tau), 1e-8)
  sigma_g <- sqrt(colSums(tau * res2) / den)

  stopifnot(!is.null(control$sigma_floor))
  sigma_g <- pmax(sigma_g, control$sigma_floor)

  em_state$beta_g <- beta_g
  em_state$eta <- mu
  em_state$sigma_g <- sigma_g
}

# Poisson M-step, no common predictors: one weighted glm.fit per component
# (weights = tau, floored at control$weight_floor). Also used for "zip".
m_step_poisson <- function(dat, em_state, control) {
  tau <- em_state$tau
  G <- em_state$G

  X <- dat$X_het
  y <- dat$y
  p <- dat$p_het

  beta_g <- matrix(NA_real_, nrow = G, ncol = p)

  for (g in seq_len(G)) {
    wg <- pmax(tau[, g], control$weight_floor)

    fit_g <- suppressWarnings(stats::glm.fit(
      x = X,
      y = y,
      weights = wg,
      family = stats::poisson()
    ))
    if (!isTRUE(fit_g$converged)){
      warning(
        "Iteration could not converge during fitting.",
        call. = FALSE
      )
    }

    coef_g <- fit_g$coefficients
    coef_g[is.na(coef_g)] <- 0

    beta_g[g, ] <- coef_g
  }

  eta <- linear_predictor_matrix(
    A = dat$X_het,
    B = dat$X_com,
    beta_g = beta_g,
    beta = NULL
  )

  em_state$beta_g <- beta_g
  em_state$eta <- eta
}

# Binomial M-step, no common predictors: one weighted glm.fit per component
# on the (successes, failures) response, using dat$binomial_size for trials.
m_step_binomial <- function(dat, em_state, control) {
  tau <- em_state$tau
  G <- em_state$G

  X <- dat$X_het
  y <- dat$y
  p <- dat$p_het

  binomial_size <- dat$binomial_size
  if (is.null(binomial_size)) {
    binomial_size <- rep(1, length(y))
  }
  binomial_size <- as.numeric(binomial_size)

  y_bin <- cbind(y, binomial_size - y)

  beta_g <- matrix(NA_real_, nrow = G, ncol = p)

  for (g in seq_len(G)) {
    wg <- pmax(tau[, g], control$weight_floor)

    fit_g <- suppressWarnings(stats::glm.fit(
      x = X,
      y = y_bin,
      weights = wg,
      family = stats::binomial()
    ))
    # This will fire multiple times with a 0/1 response. 
    # A possibility is to control how many with rlang::warn 
    if (!isTRUE(fit_g$converged)){
      warning(
        "Iteration could not converge during fitting.",
        call. = FALSE
      )
    }

    coef_g <- fit_g$coefficients
    coef_g[is.na(coef_g)] <- 0

    beta_g[g, ] <- coef_g
  }

  eta <- linear_predictor_matrix(
    A = dat$X_het,
    B = dat$X_com,
    beta_g = beta_g,
    beta = NULL
  )

  em_state$beta_g <- beta_g
  em_state$eta <- eta
}

# Gaussian M-step using structured QR

m_step_sqr_gaussian <- function(dat, em_state, control) {
  tau <- em_state$tau
  sigma_g <- em_state$sigma_g
  # Updated beta, beta_g and QR parts
  if (!is.null(sigma_g)) {
    stopifnot(length(sigma_g) == ncol(tau))

    sigma_g <- pmax(sigma_g, 1e-16)

    w_beta <- sweep(
      tau,
      MARGIN = 2,
      STATS = sigma_g^2,
      FUN = "/"
    )
  } else {
    w_beta <- tau
  }

  out <- wls_sqr(
    dat = dat,
    z = dat$y,
    weights = w_beta,
    control = control
  )
  beta <- out$beta
  beta_g <- out$beta_g

  # TODO: look at this line and think of alternatives
  mu <- linear_predictor_matrix(
    A = dat$X_het,
    B = dat$X_com,
    beta_g = beta_g,
    beta = beta
  )

  res2 <- (dat$y - mu)^2
  # Would it be faster to do these with matrix form instead of col sums?
  # TODO: test with scripts (maybe do this when in C++ land I imagine)

  den <- pmax(colSums(tau), 1e-8)
  sigma_g <- sqrt(colSums(tau * res2) / den)

  stopifnot(!is.null(control$sigma_floor))
  sigma_g <- pmax(sigma_g, control$sigma_floor)

  em_state$beta_g <- beta_g
  em_state$beta <- beta
  em_state$eta <- mu
  em_state$sigma_g <- sigma_g
}

# Poisson M-step using IRWLS and structured QR

m_step_sqr_poisson <- function(dat, em_state, control) {
  irwls_fmr(dat, em_state, control, family = "poisson")
}

# Binomial M-step using IRWLS and structured QR

m_step_sqr_binomial <- function(dat, em_state, control) {
  irwls_fmr(dat, em_state, control, "binomial")
}
