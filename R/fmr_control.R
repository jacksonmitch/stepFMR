#' Create a Control Object for stepFMR and fitFMR Models
#'
#' Creates and validates a control object of tuning parameters for
#' \code{\link{stepFMR}} and \code{\link{fitFMR}}. Most users only need
#' \code{alpha} and \code{max_iter}/\code{tol} (EM
#' convergence). The remaining parameters govern initialization strategies,
#' IRWLS loops, and numerical safeguards. They are available for tuning but
#' the defaults are reasonable for most fits. \code{fitFMR()} ignores
#' \code{alpha}, \code{parallel}, and \code{verbose},
#' which only apply to the stepwise search in \code{stepFMR()}.
#'
#' @param alpha A numeric value between 0 and 1 specifying the significance
#'  level of the hypothesis test. Default is 0.05.
#' @param verbose A logical controlling if a text progress bar is displayed
#'  during the procedure. By default is TRUE if the session is
#'  interactive, and FALSE otherwise.
#' @param parallel A logical flag for enabling user-friendly parallel
#'  computation. If \code{TRUE}, a general \code{future::multisession} plan
#'  overrides the current one, which is restored on exit. A \code{future} can be
#'  specified manually while setting this flag as \code{FALSE} for
#'  hardware-specific customization. Default is \code{FALSE}.
#' @param n_init A positive integer specifying the target number of
#'  initializations. A fixed set of initializations is used, with random starts
#'  added when \code{n_init} exceeds its size. Default is 10.
#' @param n_best_init A positive integer specifying the total number of
#'  initializations to carry into full convergence in \code{stepFMR()}.
#'  Default is 2.
#' @param init_burnin A positive integer specifying the number of EM
#'  iterations run for each candidate initialization during the burn-in
#'  stage used to pick the best starting point. Default is 10.
#' @param kmeans_starts A positive integer specifying the number of random
#'  starts to use within the K-means algorithm itself. Default is 20.
#' @param sigma_floor An optional non-negative numeric value establishing a
#'  lower bound for the estimated component standard deviations in the
#'  Gaussian family, to prevent a component from collapsing onto a single
#'  point. If \code{NULL} (the default), it is set to 5% of the standard
#'  deviation of the response. Ignored for the Poisson and binomial families.
#' @param max_iter A positive integer specifying the maximum number of
#'  iterations allowed for the algorithm. Default is 300.
#' @param tol A positive numeric value (at least 1e-10) specifying the
#'  relative change in log-likelihood below which the EM algorithm is
#'  considered converged. Default is 1e-6.
#' @param irwls_max_iter A positive integer specifying the maximum number of
#'  IRWLS iterations used for the Poisson and binomial distribution
#'  families. Default is 50.
#' @param irwls_tol A positive numeric value specifying the convergence
#'  tolerance for the IRWLS iterations used for the Poisson and binomial
#'  distribution families. Iteration stops when the largest absolute
#'  change in coefficients falls below this value. Default is 1e-6.
#'
#' @return A structured list of class \code{"FMRControl"} containing all
#'   validated control arguments and fixed internal defaults.
#'
#' @export
#'
#' @examples
#' # Generate default control options
#' ctrl <- fmr_control()
#'
#' # Customize specific hyper-parameters
#' custom_ctrl <- fmr_control(max_iter = 500, verbose = TRUE)
fmr_control <- function(
  alpha = 0.05,
  verbose = interactive(),
  parallel = FALSE,
  n_init = 10,
  n_best_init = 2,
  init_burnin = 10,
  kmeans_starts = 20,
  sigma_floor = NULL,
  max_iter = 300,
  tol = 1e-6,
  irwls_max_iter = 50,
  irwls_tol = 1e-6
) {
  # User Input checks

  collection <- checkmate::makeAssertCollection()

  checkmate::assert_number(alpha, lower = 0, upper = 1, add = collection)
  checkmate::assert_flag(verbose, add = collection)
  checkmate::assert_flag(parallel, add = collection)
  checkmate::assert_int(n_init, lower = 1, add = collection)
  checkmate::assert_int(n_best_init, lower = 1, add = collection)
  if (checkmate::test_int(n_init, lower = 1) &&
      checkmate::test_int(n_best_init, lower = 1) &&
      n_best_init > n_init) {
    collection$push(sprintf(
      "n_best_init (number of initializations carried into full EM) cannot
      exceed n_init. (received %s and %s)", n_best_init, n_init
    ))
  }
  checkmate::assert_int(init_burnin, lower = 1, add = collection)
  checkmate::assert_int(kmeans_starts, lower = 1, add = collection)
  checkmate::assert_number(sigma_floor,
    lower = 0,
    add = collection, null.ok = TRUE
  )
  checkmate::assert_int(max_iter, lower = 1, add = collection)
  checkmate::assert_number(tol, lower = 1e-10, add = collection)
  checkmate::assert_int(irwls_max_iter, lower = 1, add = collection)
  checkmate::assert_number(irwls_tol, lower = 0, add = collection)
  if (checkmate::test_number(irwls_tol) && irwls_tol == 0) {
    collection$push("irwls_tol must be greater than 0.")
  }

  if (!collection$isEmpty()) {
    err_messages <- paste0("- ", collection$getMessages(), collapse = "\n")
    stop(err_messages, call. = FALSE)
  }

  control <- list(
    alpha = alpha,
    verbose = verbose,
    parallel = parallel,
    n_init = as.integer(n_init),
    n_best_init = as.integer(n_best_init),
    init_burnin = as.integer(init_burnin),
    kmeans_starts = as.integer(kmeans_starts),
    sigma_floor = sigma_floor,
    max_iter = as.integer(max_iter),
    tol = tol,
    irwls_max_iter = as.integer(irwls_max_iter),
    irwls_tol = irwls_tol,
    init_eps = 1e-6,
    weight_floor = 1e-10
  )
  class(control) <- "FMRControl"
  control
}
