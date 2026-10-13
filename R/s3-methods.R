# ---- select_variables -------------------------------------------------------

#' Print a Variable Selection Result
#'
#' Prints a short summary of a stepwise variable-selection run: the search
#' direction, significance level, the final set of selected predictors, and
#' the resulting formula.
#'
#' @param x An object of class \code{select_variables}. This is produced
#'   internally by \code{select_variables()}, and surfaces to users as the
#'   \code{variable_selection} component of the list returned by
#'   \code{\link{stepFMR}} whenever \code{procedure} includes \code{"variables"}
#'   (see \code{stepFMR}'s \strong{Value} section).
#' @param ... Currently unused.
#'
#' @return \code{x}, invisibly.
#' @seealso \code{\link{stepFMR}}, \code{\link{print_steps.select_variables}}
#' @export
#'
#' @examples
#' betas <- rbind(
#'   c(-4, -3, 1.5, 2, 0),
#'   c(0, 1, 1.5, -1, 0),
#'   c(4, 3, 1.5, 0.5, 0)
#' )
#' pi <- c(0.2, 0.3, 0.5)
#' sigma <- c(1, 0.5, 2)
#' sim <- simulate_fmr(n = 500, betas = betas, pi = pi, sigma = sigma, seed = 1)
#' fit <- stepFMR(sim$formula, sim$data, G_max = 3, procedure = "variables")
#' # procedure = "variables" is what populates fit$variable_selection below
#' print(fit$variable_selection)
print.select_variables <- function(x, ...) {
  cat("Variable selection\n")
  cat("  Direction: ", x$direction,
    "  |  Alpha: ", x$alpha, "\n",
    sep = ""
  )
  cat("  Selected:      ", format_none(x$selected), "\n", sep = "")
  cat("  Not selected:  ", 
    format_none(setdiff(x$predictors, x$selected)), 
    "\n", sep = ""
  )
  cat("  Final formula: ", format_formula(x$final_formula), "\n", sep = "")
  invisible(x)
}

# ---- determine_effects ------------------------------------------------------

#' Print an Effect-Type Determination Result
#'
#' Prints a short summary of a stepwise effect-type determination run: the
#' search direction, significance level, which predictors ended up
#' heterogeneous vs. homogeneous, and the resulting formula.
#'
#' @param x An object of class \code{determine_effects}. This is produced
#'   internally by \code{determine_effects()}, and surfaces to users as the
#'   \code{effect_determination} component of the list returned by
#'   \code{\link{stepFMR}} whenever \code{procedure} includes \code{"effects"}
#'   (see \code{stepFMR}'s \strong{Value} section).
#' @param ... Currently unused.
#'
#' @return \code{x}, invisibly.
#' @seealso \code{\link{stepFMR}}, \code{\link{print_steps.determine_effects}}
#' @export
#'
#' @examples
#' pi <- c(0.2, 0.3, 0.5)
#' sigma <- c(1, 0.5, 2)
#' betas <- rbind(
#'   c(-4, -3, 1.5, 2),
#'   c(0, 1, 1.5, -1),
#'   c(4, 3, 1.5, 0.5)
#' )
#' sim <- simulate_fmr(n = 500, betas = betas, pi = pi, sigma = sigma, seed = 1)
#' fit <- stepFMR(sim$formula, sim$data, G_max = 3, procedure = "effects")
#' # procedure = "effects" is what populates fit$effect_determination below
#' print(fit$effect_determination)
print.determine_effects <- function(x, ...) {
  cat("Effect-type determination\n")
  cat("  Direction: ", x$direction,
    "  |  Alpha: ", x$alpha, "\n",
    sep = ""
  )
  cat("  Heterogeneous: ", format_none(x$heterogeneous), "\n", sep = "")
  cat("  Homogeneous:   ", format_none(x$homogeneous), "\n", sep = "")
  cat("  Final formula: ", format_formula(x$final_formula), "\n", sep = "")
  invisible(x)
}

# ---- print_steps ------------------------------------------------------------

#' Print the step-by-step selection table
#'
#' Generic function that prints, step by step, the candidates considered and
#' chosen at each iteration of a stepwise search. Dispatches on the class of
#' \code{x} (\code{select_variables} or \code{determine_effects}).
#'
#' @param x An object of class \code{determine_effects} or
#'   \code{select_variables}. These are exactly the \code{effect_determination}
#'   and \code{variable_selection} components of the list returned by
#'   \code{\link{stepFMR}} (present only when \code{procedure} included
#'   \code{"effects"} or \code{"variables"}, respectively).
#' @param ... Currently unused.
#'
#' @return \code{x}, invisibly.
#' @seealso \code{\link{stepFMR}}, whose \strong{Value} section documents
#'   where \code{variable_selection} and \code{effect_determination} come
#'   from.
#' @export
#'
#' @examples
#' pi <- c(0.2, 0.3, 0.5)
#' sigma <- c(1, 0.5, 2)
#' betas <- rbind(
#'   c(-4, -3, 1.5, 2, 0),
#'   c(0, 1, 1.5, -1, 0),
#'   c(4, 3, 1.5, 0.5, 0)
#' )
#' sim <- simulate_fmr(n = 500, betas = betas, pi = pi, sigma = sigma, seed = 1)
#' fit <- stepFMR(sim$formula, sim$data, G_max = 3, procedure = "variables")
#' # procedure = "variables" is what populates fit$variable_selection below
#' print_steps(fit$variable_selection)
print_steps <- function(x, ...) UseMethod("print_steps")

#' Print the Step-by-Step Effect-Type Determination Table
#'
#' For each step of the search, prints the current heterogeneous/homogeneous
#' state, every candidate's weighted p-value, and which candidate (if any)
#' was chosen.
#'
#' @param x An object of class \code{determine_effects}. This is produced
#'   internally by \code{determine_effects()}, and surfaces to users as the
#'   \code{effect_determination} component of the list returned by
#'   \code{\link{stepFMR}} whenever \code{procedure} includes \code{"effects"}
#'   (see \code{stepFMR}'s \strong{Value} section).
#' @param ... Currently unused.
#'
#' @return \code{x}, invisibly.
#' @seealso \code{\link{stepFMR}}
#' @export
#'
#' @examples
#' pi <- c(0.2, 0.3, 0.5)
#' sigma <- c(1, 0.5, 2)
#' betas <- rbind(
#'   c(-4, -3, 1.5, 2),
#'   c(0, 1, 1.5, -1),
#'   c(4, 3, 1.5, 0.5)
#' )
#' sim <- simulate_fmr(n = 500, betas = betas, pi = pi, sigma = sigma, seed = 1)
#' fit <- stepFMR(sim$formula, sim$data, G_max = 3, procedure = "effects")
#' # procedure = "effects" is what populates fit$effect_determination below
#' print_steps(fit$effect_determination)
print_steps.determine_effects <- function(x, ...) {
  label_state <- function(het, common) {
    if (length(het) == 0) {
      return("all homogeneous")
    }
    if (length(common) == 0) {
      return("all heterogeneous")
    }
    paste0(
      "heterogeneous: ", paste(het, collapse = ", "),
      "  |  homogeneous: ", paste(common, collapse = ", ")
    )
  }
  format_steps(
    x$steps, union(x$heterogeneous, x$homogeneous), x$direction, label_state
  )
  invisible(x)
}

#' Print the Step-by-Step Variable Selection Table
#'
#' For each step of the search, prints the current included/excluded state,
#' every candidate predictor's weighted p-value, and which candidate (if
#' any) was chosen.
#'
#' @param x An object of class \code{select_variables}. This is produced
#'   internally by \code{select_variables()}, and surfaces to users as the
#'   \code{variable_selection} component of the list returned by
#'   \code{\link{stepFMR}} whenever \code{procedure} includes \code{"variables"}
#'   (see \code{stepFMR}'s \strong{Value} section).
#' @param ... Currently unused.
#'
#' @return \code{x}, invisibly.
#' @seealso \code{\link{stepFMR}}
#' @export
#'
#' @examples
#' pi <- c(0.2, 0.3, 0.5)
#' sigma <- c(1, 0.5, 2)
#' betas <- rbind(
#'   c(-4, -3, 1.5, 2, 0),
#'   c(0, 1, 1.5, -1, 0),
#'   c(4, 3, 1.5, 0.5, 0)
#' )
#' sim <- simulate_fmr(n = 500, betas = betas, pi = pi, sigma = sigma, seed = 1)
#' fit <- stepFMR(sim$formula, sim$data, G_max = 3, procedure = "variables")
#' # procedure = "variables" is what populates fit$variable_selection below
#' print_steps(fit$variable_selection)
print_steps.select_variables <- function(x, ...) {
  label_state <- function(included, excluded) {
    if (length(included) == 0) {
      return("none included")
    }
    if (length(excluded) == 0) {
      return("all included")
    }
    paste0(
      "included: ", paste(included, collapse = ", "),
      "  |  excluded: ", paste(excluded, collapse = ", ")
    )
  }
  format_steps(x$steps, x$predictors, x$direction, label_state)
  invisible(x)
}

format_steps <- function(steps, all_predictors, direction, label_state) {
  in_set <- if (direction == "forward") character(0) else all_predictors

  for (s in steps) {
    out_set <- setdiff(all_predictors, in_set)
    state <- label_state(in_set, out_set)

    cat(sprintf("Step %d  (%s)\n", s$step, state))
    cat(sprintf("  %-20s  %-12s\n", "Candidate", "Weighted p-value"))
    cat(sprintf("  %-20s  %-12s\n", "---------", "-------"))

    p_fmt <- formatC(s$p0, format = "e", digits = 2)
    for (nm in names(s$p0)) {
      marker <- if (!is.na(s$chosen) && nm == s$chosen) " *" else "  "
      cat(sprintf("  %-20s  %-10s%s\n", nm, p_fmt[[nm]], marker))
    }

    if (!is.na(s$chosen)) {
      cat(sprintf("  Chosen: %s\n", s$chosen))
      in_set <- if (direction == "forward") {
        c(in_set, s$chosen)
      } else {
        setdiff(in_set, s$chosen)
      }
    } else {
      cat("  No eligible candidate -- stopping.\n")
    }
    cat("\n")
  }
}

# ---- fit_fmr ----------------------------------------------------------------

#' Print a Fitted Finite Mixture Regression Model
#'
#' Prints the number of components \code{G} and, for each parameter block
#' present on the fit, the estimated mixing proportions, component standard
#' deviations, heterogeneous coefficients, and homogeneous coefficients.
#'
#' @param x An object of class \code{fit_fmr}. This is produced internally
#'   while fitting each candidate number of components, and surfaces to
#'   users as the \code{best_fit} component of the list returned by
#'   \code{\link{stepFMR}}. It is always present, regardless of
#'   \code{procedure} (see \code{stepFMR}'s \strong{Value} section).
#' @param ... Currently unused.
#'
#' @return \code{x}, invisibly.
#' @seealso \code{\link{stepFMR}}
#' @export
#'
#' @examples
#' pi <- c(0.2, 0.3, 0.5)
#' sigma <- c(1, 0.5, 2)
#' betas <- rbind(
#'   c(-4, -3, 1.5, 2, 0),
#'   c(0, 1, 1.5, -1, 0),
#'   c(4, 3, 1.5, 0.5, 0)
#' )
#' sim <- simulate_fmr(n = 500, betas = betas, pi = pi, sigma = sigma, seed = 1)
#' fit <- stepFMR(sim$formula, sim$data, G_max = 3)
#' # fit$best_fit is always populated, at the BIC-optimal G
#' print(fit$best_fit)
print.fit_fmr <- function(x, ...) {
  cat("Best fit (G: ", x$G, ")\n", sep = "")

  if (!is.null(x$parameter_values[["pi_g"]])) {
    cat("\nMixing proportions:\n")
    print(round(x$parameter_values[["pi_g"]], 4))
  }

  if (!is.null(x$parameter_values[["sigma_g"]])) {
    cat("\nComponent standard deviations:\n")
    print(round(x$parameter_values[["sigma_g"]], 4))
  }

  if (!is.null(x$parameter_values[["beta_g"]])) {
    cat("\nHeterogeneous coefficients:\n")
    print(round(x$parameter_values[["beta_g"]], 4))
  }

  beta_present <- !is.null(x$parameter_values[["beta"]]) &&
    length(x$parameter_values[["beta"]]) > 0L
  if (beta_present) {
    cat("\nHomogeneous coefficients:\n")
    print(round(x$parameter_values[["beta"]], 4))
  }

  invisible(x)
}

# ---- stepFMR -----------------------------------------------------------------

#' Print a stepFMR Result
#'
#' Prints a top-level summary of a \code{\link{stepFMR}} fit: the family, the
#' candidate \code{G} values, and which procedure(s) were run, followed by the
#' variable-selection result (if run), the effect-determination result (if
#' run), and the best fit.
#'
#' @param x An object of class \code{stepFMR}, as returned by
#'   \code{\link{stepFMR}}.
#' @param ... Currently unused.
#'
#' @return \code{x}, invisibly.
#' @export
#'
#' @examples
#' pi <- c(0.2, 0.3, 0.5)
#' sigma <- c(1, 0.5, 2)
#' betas <- rbind(
#'   c(-4, -3, 1.5, 2, 0),
#'   c(0, 1, 1.5, -1, 0),
#'   c(4, 3, 1.5, 0.5, 0)
#' )
#' sim <- simulate_fmr(n = 500, betas = betas, pi = pi, sigma = sigma, seed = 1)
#' fit <- stepFMR(sim$formula, sim$data, G_max = 3)
#' print(fit)
print.stepFMR <- function(x, ...) {
  cat("stepFMR | family: ", x$family,
    " | G: ", paste(x$G_values, collapse = ","),
    " | procedure: ", x$procedure, "\n\n",
    sep = ""
  )

  if (!is.null(x$variable_selection)) {
    print(x$variable_selection)
    cat("\n")
  }

  if (!is.null(x$effect_determination)) {
    print(x$effect_determination)
    cat("\n")
  }

  if (is.null(x$effect_determination) &&
        !is.null(x$variable_selection) &&
        length(x$variable_selection$selected) == 0L) {
    cat("Effect-type determination skipped: no predictors selected.\n\n")
  }

  if (!is.null(x$best_fit)) {
    print(x$best_fit)
  }

  invisible(x)
}

#' Summarize a stepFMR Result
#'
#' Collects the call, formula, family, procedure, candidate \code{G} values, best
#' fit, and (if run) the variable-selection and effect-determination results
#' and their step-by-step tables into a single \code{summary.stepFMR} object
#' for printing.
#'
#' @param object An object of class \code{stepFMR}, as returned by
#'   \code{\link{stepFMR}}.
#' @param ... Currently unused.
#'
#' @return An object of class \code{summary.stepFMR}: a list containing the
#'   call, formula, family, procedure, candidate \code{G} values, the best fit,
#'   the variable-selection and effect-determination results (if run), and
#'   their step-by-step tables.
#' @export
#'
#' @examples
#' pi <- c(0.2, 0.3, 0.5)
#' sigma <- c(1, 0.5, 2)
#' betas <- rbind(
#'   c(-4, -3, 1.5, 2, 0),
#'   c(0, 1, 1.5, -1, 0),
#'   c(4, 3, 1.5, 0.5, 0)
#' )
#' sim <- simulate_fmr(n = 500, betas = betas, pi = pi, sigma = sigma, seed = 1)
#' fit <- stepFMR(sim$formula, sim$data, G_max = 3)
#' summary <- summary(fit)
summary.stepFMR <- function(object, ...) {
  out <- list(
    call = object$call,
    formula = object$formula,
    family = object$family,
    procedure = object$procedure,
    G_values = object$G_values,
    best_fit = object$best_fit,
    selected_G = object$best_fit$G,
    variable_selection = object$variable_selection,
    variable_steps = object$variable_selection$steps,
    effect_determination = object$effect_determination,
    effect_steps = object$effect_determination$steps
  )
  class(out) <- "summary.stepFMR"
  out
}

#' Print a stepFMR Summary
#'
#' Prints the top-level summary produced by \code{\link{summary.stepFMR}},
#' including the step-by-step tables for the variable-selection and
#' effect-determination searches (if run), and the best fit.
#'
#' @param x An object of class \code{summary.stepFMR}, as returned by
#'   \code{\link{summary.stepFMR}}.
#' @param ... Currently unused.
#'
#' @return \code{x}, invisibly.
#' @export
#'
#' @examples
#' pi <- c(0.2, 0.3, 0.5)
#' sigma <- c(1, 0.5, 2)
#' betas <- rbind(
#'   c(-4, -3, 1.5, 2, 0),
#'   c(0, 1, 1.5, -1, 0),
#'   c(4, 3, 1.5, 0.5, 0)
#' )
#' sim <- simulate_fmr(n = 500, betas = betas, pi = pi, sigma = sigma, seed = 1)
#' fit <- stepFMR(sim$formula, sim$data, G_max = 3)
#' summary(fit)
#' print(summary(fit))
print.summary.stepFMR <- function(x, ...) {
  cat("stepFMR summary | family: ", x$family,
    " | G: ", paste(x$G_values, collapse = ","),
    " | procedure: ", x$procedure, "\n\n",
    sep = ""
  )

  if (!is.null(x$variable_selection)) {
    print(x$variable_selection)
    cat("\n")
    print_steps(x$variable_selection)
  }

  if (!is.null(x$effect_determination)) {
    print(x$effect_determination)
    cat("\n")
    print_steps(x$effect_determination)
  }

  if (is.null(x$effect_determination) &&
        !is.null(x$variable_selection) &&
        length(x$variable_selection$selected) == 0L) {
    cat("Effect-type determination skipped: no predictors selected.\n\n")
  }

  if (!is.null(x$best_fit)) {
    print(x$best_fit)
  }

  invisible(x)
}

# ---- fitFMR ------------------------------------------------------------------

#' Print a fitFMR Result
#'
#' Prints the family, formula, and which predictors are heterogeneous or
#' homogeneous, a table of log-likelihood, number of parameters, and BIC for
#' each fitted number of components (the selected one is marked with
#' \code{*}), and the estimates of the selected fit.
#'
#' @param x An object of class \code{fitFMR}, as returned by
#'   \code{\link{fitFMR}}.
#' @param ... Currently unused.
#'
#' @return \code{x}, invisibly.
#' @seealso \code{\link{fitFMR}}, \code{\link{summary.fitFMR}}
#' @export
#'
#' @examples
#' sim <- fmr_gaussian_example
#' fit <- fitFMR(sim$formula, sim$data, G = 2:3, common = "x2")
#' print(fit)
print.fitFMR <- function(x, ...) {
  cat("fitFMR | family: ", x$family,
    " | G: ", paste(x$G_values, collapse = ","), "\n",
    sep = ""
  )
  cat("  Formula:       ", format_formula(x$formula), "\n", sep = "")
  cat("  Heterogeneous: ", format_none(x$heterogeneous), "\n", sep = "")
  cat("  Homogeneous:   ", format_none(x$common), "\n\n", sep = "")

  tab <- x$table[, c("G", "loglik", "num_parameters", "bic")]
  tab$loglik <- round(tab$loglik, 2)
  tab$bic <- round(tab$bic, 2)
  tab[[" "]] <- ifelse(tab$G == x$selected_G, "*", "")
  print(tab, row.names = FALSE)
  cat("\n")

  print(x$best_fit)
  invisible(x)
}

#' Summarize a fitFMR Result
#'
#' Collects the call, formula, family, effect types, the table of fits for
#' each number of components, the selected fit, and its initialization
#' diagnostics into a \code{summary.fitFMR} object for printing.
#'
#' @param object An object of class \code{fitFMR}, as returned by
#'   \code{\link{fitFMR}}.
#' @param ... Currently unused.
#'
#' @return An object of class \code{summary.fitFMR}: a list containing the
#'   call, formula, family, candidate \code{G} values, common and
#'   heterogeneous predictors, the table of fits, the selected \code{G}, the
#'   selected fit, and its initialization diagnostics.
#' @seealso \code{\link{fitFMR}}, \code{\link{print.fitFMR}}
#' @export
#'
#' @examples
#' sim <- fmr_gaussian_example
#' fit <- fitFMR(sim$formula, sim$data, G = 2:3, common = "x2")
#' summary(fit)
summary.fitFMR <- function(object, ...) {
  out <- list(
    call = object$call,
    formula = object$formula,
    family = object$family,
    G_values = object$G_values,
    common = object$common,
    heterogeneous = object$heterogeneous,
    table = object$table,
    selected_G = object$selected_G,
    best_fit = object$best_fit,
    initialization = object$best_fit$initialization
  )
  class(out) <- "summary.fitFMR"
  out
}

#' Print a fitFMR Summary
#'
#' Prints the summary produced by \code{\link{summary.fitFMR}}: the model
#' structure, the full table of fits (including iterations and convergence),
#' a one-line initialization summary for the selected fit, and the selected
#' fit's estimates.
#'
#' @param x An object of class \code{summary.fitFMR}, as returned by
#'   \code{\link{summary.fitFMR}}.
#' @param ... Currently unused.
#'
#' @return \code{x}, invisibly.
#' @seealso \code{\link{fitFMR}}
#' @export
#'
#' @examples
#' sim <- fmr_gaussian_example
#' fit <- fitFMR(sim$formula, sim$data, G = 2:3, common = "x2")
#' print(summary(fit))
print.summary.fitFMR <- function(x, ...) {
  cat("fitFMR summary | family: ", x$family,
    " | G: ", paste(x$G_values, collapse = ","), "\n",
    sep = ""
  )
  cat("  Formula:       ", format_formula(x$formula), "\n", sep = "")
  cat("  Heterogeneous: ", format_none(x$heterogeneous), "\n", sep = "")
  cat("  Homogeneous:   ", format_none(x$common), "\n\n", sep = "")

  tab <- x$table
  tab$loglik <- round(tab$loglik, 2)
  tab$bic <- round(tab$bic, 2)
  tab[[" "]] <- ifelse(tab$G == x$selected_G, "*", "")
  print(tab, row.names = FALSE)
  cat("\n")

  init <- x$initialization
  cat(sprintf(
    "Initialization (G = %d): %d starts tried, %d kept; best start: %s\n",
    x$selected_G, init$n_starts, init$n_kept, x$best_fit$best_init_name
  ))
  if (any(!x$table$converged)) {
    cat("Note: EM did not converge for G = ",
      paste(x$table$G[!x$table$converged], collapse = ", "), "\n",
      sep = ""
    )
  }
  cat("\n")

  print(x$best_fit)
  invisible(x)
}
