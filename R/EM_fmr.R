# EM algorithm for finite mixture regression

#' Run the EM Algorithm to Convergence
#'
#' Alternates the M-step (\code{m_step()}) and E-step (\code{e_step_fmr()})
#' until the change in log-likelihood falls below \code{control$tol} or
#' \code{control$max_iter} is reached. \code{em_state} is mutated in place
#' over the course of the iterations.
#'
#' @param prepared_data A \code{FMRData} object (from \code{prepare_data()}).
#' @param G An integer number of mixture components.
#' @param em_state An \code{EmState} object holding the initial (and,
#'   afterward, final) EM state.
#' @param family A character string specifying the error distribution:
#'   \code{"gaussian"}, \code{"poisson"}, or \code{"binomial"}.
#' @param control A \code{FMRControl} object (from \code{fmr_control()}).
#'
#' @return A list with elements: \code{em_state} (the final, mutated
#'   \code{EmState}), \code{loglik} (final log-likelihood),
#'   \code{loglik_trace} (log-likelihood at each iteration),
#'   \code{iterations} (number of iterations run), \code{converged}
#'   (whether the tolerance criterion was met), and
#'   \code{irwls_iterations}/\code{irwls_converged} (per-iteration IRWLS
#'   diagnostics; \code{NULL} for the Gaussian family).
#' @noRd
em_fmr <- function(em_state, prepared_data, G, family, control) {
  stopifnot(inherits(em_state, "EmState"))

  loglik_trace <- numeric(control$max_iter)
  irwls_iterations <- integer(control$max_iter)
  irwls_converged <- logical(control$max_iter)
  converged <- FALSE
  loglik_old <- -Inf

  for (iter in seq_len(control$max_iter)) {
    m_step(prepared_data, em_state, family, control)

    irwls_iterations[iter] <- em_state$irwls_iterations
    irwls_converged[iter] <- em_state$irwls_converged

    e_step_fmr(prepared_data, em_state, family)

    loglik_new <- em_state$loglik
    loglik_trace[iter] <- loglik_new

    tol <- control$tol
    if (is.finite(loglik_old) && tol > 0) {
      ll_diff <- abs(loglik_new - loglik_old)
      ll_scale <- 1 + abs(loglik_old)
      if (ll_diff < tol * ll_scale) {
        converged <- TRUE
        break
      }
    }

    loglik_old <- loglik_new
  }

  loglik_trace <- loglik_trace[seq_len(iter)]
  irwls_iterations <- irwls_iterations[seq_len(iter)]
  irwls_converged <- irwls_converged[seq_len(iter)]

  list(
    em_state = em_state,
    loglik = loglik_trace[length(loglik_trace)],
    loglik_trace = loglik_trace,
    iterations = iter,
    converged = converged,
    irwls_iterations = if (family == "gaussian") NULL else irwls_iterations,
    irwls_converged = if (family == "gaussian") NULL else irwls_converged
  )
}
