#' Mutable EM Iteration State
#'
#' An R6 class holding the mutable state of one EM run: responsibilities
#' \code{tau}, mixing proportions \code{pi_g}, heterogeneous/common
#' coefficients \code{beta_g}/\code{beta}, Gaussian component standard
#' deviations \code{sigma_g}, the linear predictor \code{eta}, the current
#' \code{loglik}, and IRWLS diagnostics. Passed by reference into
#' \code{m_step()}/\code{e_step()}/\code{irwls_fmr()}, which mutate it
#' in place rather than returning a new object each iteration.
#'
#' @param tau A numeric responsibility matrix (n x G).
#' @param pi_g An optional numeric vector of mixing proportions; defaults to
#'   \code{colMeans(tau)}.
#' @param beta_g An optional numeric matrix of heterogeneous coefficients
#'   (G x p).
#' @param beta An optional numeric vector of common/homogeneous
#'   coefficients.
#' @param sigma_g An optional numeric vector of component standard
#'   deviations (Gaussian family only).
#' @param eta An optional numeric matrix of linear predictor values.
#' @param loglik A numeric log-likelihood value. Defaults to \code{-Inf}.
#' @param irwls_iterations An integer count of IRWLS iterations run in the
#'   most recent M-step (non-Gaussian families only).
#' @param irwls_converged A logical flag for whether IRWLS converged in the
#'   most recent M-step (non-Gaussian families only).
#'
#' @return A new \code{EmState} object.
#' @noRd
EmState <- R6::R6Class("EmState",
  private = list(
    .tau              = NULL,
    .pi_g             = NULL,
    .beta_g           = NULL,
    .beta             = NULL,
    .sigma_g          = NULL,
    .eta              = NULL,
    .loglik           = NULL,
    .irwls_iterations = NULL,
    .irwls_converged  = NULL
  ),
  active = list(
    tau = function(value) {
      if (missing(value)) {
        return(private$.tau)
      }
      stopifnot(is.matrix(value))
      private$.tau <- value
    },
    G = function(value) {
      if (missing(value)) {
        return(ncol(private$.tau))
      }
      stop("G is read-only; it is derived from tau.")
    },
    pi_g = function(value) {
      if (missing(value)) {
        return(private$.pi_g)
      }
      stopifnot(is.numeric(value))
      private$.pi_g <- value
    },
    # These parameters can be set/start as NULL
    beta_g = function(value) {
      if (missing(value)) {
        return(private$.beta_g)
      }
      stopifnot(is.null(value) || is.matrix(value))
      private$.beta_g <- value
    },
    beta = function(value) {
      if (missing(value)) {
        return(private$.beta)
      }
      stopifnot(is.null(value) || is.numeric(value))
      private$.beta <- value
    },
    sigma_g = function(value) {
      if (missing(value)) {
        return(private$.sigma_g)
      }
      stopifnot(is.null(value) || is.numeric(value))
      private$.sigma_g <- value
    },
    eta = function(value) {
      if (missing(value)) {
        return(private$.eta)
      }
      stopifnot(is.null(value) || is.numeric(value))
      private$.eta <- value
    },
    loglik = function(value) {
      if (missing(value)) {
        return(private$.loglik)
      }
      stopifnot(is.numeric(value), length(value) == 1)
      private$.loglik <- value
    },
    irwls_iterations = function(value) {
      if (missing(value)) {
        return(private$.irwls_iterations)
      }
      stopifnot(is.numeric(value), length(value) == 1)
      private$.irwls_iterations <- value
    },
    irwls_converged = function(value) {
      if (missing(value)) {
        return(private$.irwls_converged)
      }
      stopifnot(is.logical(value), length(value) == 1)
      private$.irwls_converged <- value
    }
  ),
  public = list(
    initialize = function(tau,
                          pi_g = NULL,
                          beta_g = NULL,
                          beta = NULL,
                          sigma_g = NULL,
                          eta = NULL,
                          loglik = -Inf,
                          irwls_iterations = NA_integer_,
                          irwls_converged = FALSE) {
      self$tau <- tau
      self$pi_g <- if (is.null(pi_g)) colMeans(tau) else pi_g

      self$beta_g <- beta_g
      self$beta <- beta
      self$sigma_g <- sigma_g
      self$eta <- eta
      self$loglik <- loglik

      self$irwls_iterations <- irwls_iterations
      self$irwls_converged <- irwls_converged
    },
    to_list = function(prepared_data = NULL) {
      fields <- list(
        tau = self$tau,
        pi_g = self$pi_g,
        beta_g = self$beta_g,
        beta = self$beta,
        sigma_g = self$sigma_g,
        eta = self$eta,
        loglik = self$loglik
      )

      fields <- compact(fields)

      if (!is.null(prepared_data)) {
        if (!is.null(fields[["beta_g"]])) {
          rownames(fields[["beta_g"]]) <- paste0(
            "g",
            seq_len(nrow(fields[["beta_g"]]))
          )
          if (!is.null(colnames(prepared_data$X_het))) {
            colnames(fields[["beta_g"]]) <- colnames(prepared_data$X_het)
          }
        }
        if (!is.null(fields[["beta"]]) && length(fields[["beta"]]) > 0L) {
          fields[["beta"]] <- as.numeric(fields[["beta"]])
          if (!is.null(colnames(prepared_data$X_com))) {
            names(fields[["beta"]]) <- colnames(prepared_data$X_com)
          }
        }
        if (!is.null(fields[["sigma_g"]])) {
          names(fields[["sigma_g"]]) <- paste0(
            "g",
            seq_along(fields[["sigma_g"]])
          )
        }
        if (!is.null(fields[["pi_g"]])) {
          names(fields[["pi_g"]]) <- paste0(
            "g",
            seq_along(fields[["pi_g"]])
          )
        }
        if (!is.null(fields[["tau"]])) {
          colnames(fields[["tau"]]) <- paste0(
            "g",
            seq_len(ncol(fields[["tau"]]))
          )
        }
      }

      fields
    }
  )
)
