# Run the WEST procedure for one candidate predictor.
#
# shared_fits : output of fit_across_G for the shared side (null for forward,
#               alt for backward).
# direction   : "forward" or "backward",determines which role shared_fits plays
# candidate_data : prepared_data for the candidate side (already built
#                           by the caller for this specific predictor).

#' Compare a Candidate Predictor Against a Shared Baseline Fit
#'
#' For one candidate predictor, fits the candidate model across all values
#' of \code{model$G_values} (warm-started from the shared fit's
#' responsibilities) and compares it to the shared/baseline fit at each
#' \code{G} via a likelihood ratio test. The per-\code{G} p-values are then
#' combined into a single BIC-weighted p-value used for the stepwise
#' eligibility decision.
#'
#' @param model A \code{WMRModel} object.
#' @param direction A character string, either \code{"forward"} or
#'   \code{"backward"}, determining whether \code{shared_fits} plays the
#'   role of the null or the alternative model in the likelihood ratio test.
#' @param shared_fits A list of \code{fit_fmr} objects (one per value in
#'   \code{model$G_values}) for the shared/baseline model.
#' @param candidate_data A \code{WMRData} object (from \code{prepare_data()})
#'   for the candidate predictor's model specification.
#'
#' @return A list of class \code{west_procedure} with elements: \code{p0}
#'   (the BIC-weighted p-value), \code{reject} (whether \code{p0 < alpha}),
#'   \code{alpha}, \code{table} (a data frame of per-\code{G} log-likelihoods,
#'   BIC values, LRT statistics, degrees of freedom, p-values, and BIC
#'   \code{weights}), \code{shared_fits} (the baseline fits across
#'   \code{model$G_values}), and \code{candidate_fits} (fits across
#'   \code{model$G_values} for the candidate model).
#' @noRd
west_procedure <- function(
  model,
  direction,
  shared_fits,
  candidate_data
) {
  G_values <- model$G_values
  alpha <- model$control$alpha

  stopifnot(length(shared_fits) == length(G_values))

  candidate_fits <- fit_across_G(model, candidate_data,
     init_lists = make_warm_init_lists(
      shared_fits, "previous_fit", model, candidate_data))

  rows <- vector("list", length(G_values))

  for (i in seq_along(G_values)) {
    G <- G_values[[i]]

    if (direction == "forward") {
      fit_null <- shared_fits[[i]]
      fit_alt <- candidate_fits[[i]]
    } else {
      fit_null <- candidate_fits[[i]]
      fit_alt <- shared_fits[[i]]
    }

    k_null <- fit_null$num_parameters
    k_alt <- fit_alt$num_parameters

    lrt <- -2 * (fit_null$loglik - fit_alt$loglik)
    df <- k_alt - k_null

    bic_null <- fit_null$bic
    bic_alt <- fit_alt$bic

    p_value <- if (is.finite(lrt) && is.finite(df) && df > 0) {
      stats::pchisq(lrt, df = df, lower.tail = FALSE)
    } else {
      NA_real_
    }

    rows[[i]] <- data.frame(
      G = G,
      loglik_null = fit_null$loglik,
      loglik_alt = fit_alt$loglik,
      bic_null = bic_null,
      bic_alt = bic_alt,
      lrt = lrt,
      df = df,
      p_value = p_value,
      null_converged = fit_null$converged,
      alt_converged = fit_alt$converged,
      stringsAsFactors = FALSE
    )
  }

  tab <- do.call(rbind, rows)

  if (direction == "forward") {
    model_bic <- tab$bic_null
  } else {
    model_bic <- tab$bic_alt
  }
  valid <- is.finite(model_bic) & is.finite(tab$p_value)

  if (!any(valid)) {
    p0 <- NA_real_
    weights <- rep(NA_real_, nrow(tab))
    reject <- FALSE
  } else {
    bic_valid <- model_bic[valid]
    shifted <- -0.5 * (bic_valid - min(bic_valid))
    w_valid <- exp(shifted) / sum(exp(shifted))

    weights <- rep(NA_real_, nrow(tab))
    weights[valid] <- w_valid

    p0 <- sum(tab$p_value[valid] * w_valid)
    reject <- is.finite(p0) && p0 < alpha
  }

  tab$weights <- weights

  out <- list(
    p0 = p0,
    reject = reject,
    alpha = alpha,
    table = tab,
    shared_fits = shared_fits,
    candidate_fits = candidate_fits
  )

  class(out) <- "west_procedure"
  out
}
