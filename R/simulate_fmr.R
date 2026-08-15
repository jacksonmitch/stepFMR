#' Simulate Data from a Finite Mixture Regression
#'
#' Simulates covariates, group membership, and a response from a specified
#' finite mixture regression model. Column/row names on `betas` are
#' optional. If omitted, columns become `"(Intercept)", "x1", "x2", ...`
#' and rows become `"g1", "g2", ...`.
#'
#' @param n A positive integer, the number of observations to simulate.
#' @param betas A numeric matrix of true coefficients, one row per mixture
#'   component (or a list of numeric vectors, one per component, which
#'   will be row-bound into a matrix). If column names are supplied, one
#'   of them must be `"(Intercept)"`; if omitted, the first column is
#'   assumed to be the intercept and named accordingly.
#' @param pi An optional numeric vector of mixing proportions, one per
#'   component (normalized internally). Defaults to equal proportions.
#' @param sigma A numeric vector of per-component residual standard
#'   deviations (Gaussian family only). Recycled to the number of
#'   components if length 1. Default is 1.
#' @param family A character string, one of `"gaussian"`, `"poisson"`, or
#'   `"binomial"`.
#' @param size For `family = "binomial"`, the number of trials per
#'   observation (recycled if length 1). Default is 1.
#' @param seed An optional integer passed to `set.seed()` for reproducible
#'   simulation.
#'
#' @return A list with elements: `data` (a data.frame with the response
#'   `y`, the covariates, and — for `family = "binomial"` — a `size` trial
#'   count column), `formula` (ready to pass to `stepFMR()`; a
#'   `cbind(y, size - y) ~ .` form for binomial), and `true_group` (the
#'   simulated integer group labels, for use with `match_groups()`).
#' @export
#'
#' @examples
#' betas <- rbind(
#'   c(-2, -3, 1),
#'   c(2, 3, 1)
#' )
#' sim <- simulate_fmr(n = 300, betas = betas, seed = 1)
#' fit <- stepFMR(sim$formula, sim$data, G_max = 3)
#' fit_aligned <- match_groups(fit, sim$true_group)
simulate_fmr <- function(n,
                         betas,
                         pi = NULL,
                         sigma = 1,
                         family = c("gaussian", "poisson", "binomial"),
                         size = 1,
                         seed = NULL) {
  family <- match.arg(family)

  checkmate::assert_count(n, positive = TRUE)

  if (is.list(betas) && !is.matrix(betas)) {
    betas <- do.call(rbind, betas)
  }
  checkmate::assert_matrix(betas, mode = "numeric")

  if (is.null(colnames(betas))) {
    colnames(betas) <- c("(Intercept)", paste0("x", seq_len(ncol(betas) - 1)))
  } else if (!"(Intercept)" %in% colnames(betas)) {
    stop("betas has column names but none is 'Intercept'.", call. = FALSE)
  }
  if (is.null(rownames(betas))) {
    rownames(betas) <- paste0("g", seq_len(nrow(betas)))
  }

  if (!is.null(seed)) set.seed(seed)

  G <- nrow(betas)

  if (!is.null(pi)) {
    checkmate::assert_numeric(pi, len = G, lower = 0)
  } else {
    pi <- rep(1 / G, G)
  }
  checkmate::assert_numeric(sigma, lower = 0)
  sigma <- rep_len(sigma, G)
  checkmate::assert_numeric(size, lower = 1)

  z <- sample.int(G, size = n, replace = TRUE, prob = pi / sum(pi))

  pred_names <- setdiff(colnames(betas), "(Intercept)")
  X <- matrix(stats::rnorm(n * length(pred_names)),
    nrow = n,
    dimnames = list(NULL, pred_names)
  )
  Xfull <- cbind("(Intercept)" = 1, X)

  eta <- rowSums(Xfull * betas[z, colnames(Xfull), drop = FALSE])

  size_vec <- if (family == "binomial") rep_len(size, n) else NULL

  y <- switch(family,
    gaussian = stats::rnorm(n, mean = eta, sd = sigma[z]),
    poisson = stats::rpois(n, lambda = exp(eta)),
    binomial = stats::rbinom(n, size = size_vec, prob = stats::plogis(eta))
  )

  data <- cbind(data.frame(y = y), as.data.frame(X))

  if (family == "binomial") {
    data$size <- size_vec
    formula <- stats::as.formula(paste(
      "cbind(y, size - y) ~", paste(pred_names, collapse = " + ")
    ))
  } else {
    formula <- make_formula(pred_names, "y")
  }

  list(data = data, formula = formula, true_group = z)
}