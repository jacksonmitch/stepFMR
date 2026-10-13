#' Number of Free Parameters
#'
#' Counts the free parameters of a fitted mixture regression: component
#' coefficients, common coefficients, mixing proportions (\code{G - 1}),
#' and, for Gaussian, component standard deviations. For \code{"zip"} only
#' \code{G - 1} components carry coefficients.
#'
#' @param ncol_het,ncol_common Number of heterogeneous and common columns in
#'   the design.
#' @param G Number of mixture components.
#' @param family Response family.
#'
#' @return An integer count.
#' @noRd
count_params_fmr <- function(ncol_het, ncol_common, G, family) {

  k <- if (family == "zip") {
    (G - 1) * ncol_het + ncol_common + (G - 1)
  } else if (family == "gaussian") {
    G * ncol_het + ncol_common + (G - 1) + G
  } else {
    G * ncol_het + ncol_common + (G - 1)
  }

  k
}

# BIC = -2 * loglik + k * log(n). Smaller is better.
compute_bic <- function(loglik, n, k) {
  -2 * loglik + k * log(n)
}

# Pick the fit with the best criterion value (lowest BIC, or highest
# loglik) from a list of fit_fmr objects. Non-finite values are ignored.
select_best_G <- function(fits, criterion = c("bic", "loglik")) {
  criterion <- match.arg(criterion)

  values <- vapply(fits, function(fit) {
    value <- fit[[criterion]]

    if (is.null(value) || length(value) != 1L) {
      return(NA_real_)
    }

    as.numeric(value)
  }, numeric(1))

  valid <- is.finite(values)

  if (criterion == "bic") {
    best_index <- which.min(ifelse(valid, values, Inf))
  } else {
    best_index <- which.max(ifelse(valid, values, -Inf))
  }

  best_fit <- fits[[best_index]]
  best_fit
}
