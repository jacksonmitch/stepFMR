# Validated data container for mixture regression fitting.
# All structural checks happen once here so callers can trust the contents.

#' Build a Validated Design-Matrix Container
#'
#' Constructs the heterogeneous and common design matrices for a given model
#' specification from \code{model$formula}/\code{model$data}, and wraps them
#' (with the response and, for binomial models, the binomial size) in a
#' validated \code{FMRData} object.
#'
#' @param model A \code{FMRModel} object.
#' @param included A character vector of predictor names to include in the
#'   design. Defaults to all predictors in \code{model}.
#' @param common A character vector of predictor names, a subset of
#'   \code{included}, to treat as common/homogeneous (i.e. excluded from the
#'   heterogeneous design matrix). Defaults to \code{NULL} (all included
#'   predictors are heterogeneous).
#'
#' @return A \code{FMRData} object.
#' @noRd
prepare_data <- function(model, included = model$predictors, common = NULL) {
  mf <- model$mf
  y <- as.numeric(stats::model.response(mf))

  binomial_size_vec <- NULL
  if (identical(model$family, "binomial")) {
    stopifnot(".binom_size" %in% names(model$data))
    binomial_size_vec <- as.numeric(model$data$.binom_size)
  }

  het_formula <- make_formula(setdiff(included, common), response = model$response)
  X_het <- stats::model.matrix(het_formula, data = mf)

  if (length(common) > 0) {
    X_com <- stats::model.matrix(
      stats::as.formula(paste("~", paste(common, collapse = " + "), "- 1")), 
      data = mf)
  } else {
    X_com <- matrix(numeric(0), nrow = length(y), ncol = 0)
  }

  FMRData$new(y = y, X_het = X_het, X_com = X_com, binomial_size = binomial_size_vec)
}

#' Validated Design-Matrix Container
#'
#' An R6 class holding one validated design-matrix specification for
#' mixture regression fitting: the response \code{y}, the heterogeneous
#' design matrix \code{X_het}, the common design matrix \code{X_com}, and
#' (for binomial models) \code{binomial_size}. All fields are read-only
#' active bindings; structural checks (matching dimensions, no missing
#' values, \code{n > ncol(X_het)}) run once in \code{initialize()}.
#'
#' @param y A numeric response vector.
#' @param X_het A numeric matrix, the heterogeneous/component-specific
#'   design matrix.
#' @param X_com A numeric matrix, the common/homogeneous design matrix.
#' @param binomial_size An optional numeric vector of binomial trial counts.
#'
#' @return A new \code{FMRData} object.
#' @noRd
FMRData <- R6::R6Class(
  "FMRData",
  private = list(
    .y = NULL,
    .X_het = NULL,
    .X_com = NULL,
    .binomial_size = NULL
  ),
  active = list(
    y = function(v) {
      if (missing(v)) private$.y else stop("y is read-only")
    },
    X_het = function(v) {
      if (missing(v)) private$.X_het else stop("X_het is read-only")
    },
    X_com = function(v) {
      if (missing(v)) private$.X_com else stop("X_com is read-only")
    },
    binomial_size = function(v) {
      if (missing(v)) private$.binomial_size else stop("binomial_size is read-only")
    },
    n = function() nrow(private$.X_het),
    p_het = function() ncol(private$.X_het),
    p_com = function() ncol(private$.X_com)
  ),
  public = list(
    initialize = function(y, X_het, X_com, binomial_size = NULL) {
      if (!is.numeric(y)) stop("y must be a numeric vector.")
      if (!is.matrix(X_het)) stop("X_het must be a matrix.")
      if (!is.matrix(X_com)) stop("X_com must be a matrix.")

      n <- length(y)

      if (nrow(X_het) != n) {
        stop(sprintf(
          "nrow(X_het) (%d) must equal length(y) (%d).",
          nrow(X_het), n
        ))
      }
      if (nrow(X_com) != n) {
        stop(sprintf(
          "nrow(X_com) (%d) must equal length(y) (%d).",
          nrow(X_com), n
        ))
      }

      if (n <= ncol(X_het)) {
        stop(sprintf(
          "Need n > ncol(X_het) for the weighted QR M-step.
          Got n = %d, ncol(X_het) = %d.",
          n, ncol(X_het)
        ))
      }

      # na.fail in model.frame covers y and X_het; X_com is built separately
      if (anyNA(X_com)) {
        stop("X_com contains missing values.")
      }
      private$.y <- y
      private$.X_het <- X_het
      private$.X_com <- X_com
      private$.binomial_size <- binomial_size
    },
    print = function(...) {
      cat(
        "<FMRData>  n =", self$n,
        " p_het =", self$p_het,
        " p_com =", self$p_com, "\n"
      )
      invisible(self)
    }
  )
)
