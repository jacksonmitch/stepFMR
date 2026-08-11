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
#'   \code{\link{westMR}} whenever \code{procedure} includes \code{"variables"}
#'   (see \code{westMR}'s \strong{Value} section).
#' @param ... Currently unused.
#'
#' @return \code{x}, invisibly.
#' @seealso \code{\link{westMR}}, \code{\link{steps_table.select_variables}}
#' @export
#'
#' @examples
#' set.seed(1)
#' n <- 500
#' pi <- c(0.2, 0.3, 0.5)
#' z <- sample(1:3, n, replace = TRUE, prob = pi)
#' x1 <- rnorm(n)
#' x2 <- rnorm(n)
#' x3 <- rnorm(n)
#' x4 <- rnorm(n)
#' sigma <- c(1, 0.5, 2)
#' beta <- rbind(
#'   g1 = c(intercept = -4, x1 = -3, x2 = 1.5, x3 = 2, x4 = 0),
#'   g2 = c(intercept = 0, x1 = 1, x2 = 1.5, x3 = -1, x4 = 0),
#'   g3 = c(intercept = 4, x1 = 3, x2 = 1.5, x3 = 0.5, x4 = 0)
#' )
#' eta <- beta[z, "intercept"] + beta[z, "x1"] * x1 + beta[z, "x2"] * x2 +
#'   beta[z, "x3"] * x3 + beta[z, "x4"] * x4
#' y <- rnorm(n, mean = eta, sd = sigma[z])
#' dat <- data.frame(y = y, x1 = x1, x2 = x2, x3 = x3, x4 = x4)
#'
#' # x4 has no true effect and should be dropped during selection
#' fit <- westMR(
#'   y ~ x1 + x2 + x3 + x4, data = dat, G_max = 3, procedure = "variables"
#' )
#' # procedure = "variables" is what populates fit$variable_selection below
#' print(fit$variable_selection)
print.select_variables <- function(x, ...) {
  cat("Variable selection\n")
  cat("  Direction: ", x$direction,
    "  |  Alpha: ", x$alpha, "\n",
    sep = ""
  )
  cat("  Selected:      ", format_none(x$selected), "\n", sep = "")
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
#'   \code{\link{westMR}} whenever \code{procedure} includes \code{"effects"}
#'   (see \code{westMR}'s \strong{Value} section).
#' @param ... Currently unused.
#'
#' @return \code{x}, invisibly.
#' @seealso \code{\link{westMR}}, \code{\link{steps_table.determine_effects}}
#' @export
#'
#' @examples
#' set.seed(1)
#' n <- 500
#' pi <- c(0.2, 0.3, 0.5)
#' z <- sample(1:3, n, replace = TRUE, prob = pi)
#' x1 <- rnorm(n)
#' x2 <- rnorm(n)
#' x3 <- rnorm(n)
#' sigma <- c(1, 0.5, 2)
#' beta <- rbind(
#'   g1 = c(intercept = -4, x1 = -3, x2 = 1.5, x3 = 2),
#'   g2 = c(intercept = 0, x1 = 1, x2 = 1.5, x3 = -1),
#'   g3 = c(intercept = 4, x1 = 3, x2 = 1.5, x3 = 0.5)
#' )
#' eta <- beta[z, "intercept"] + beta[z, "x1"] * x1 + beta[z, "x2"] * x2 +
#'   beta[z, "x3"] * x3
#' y <- rnorm(n, mean = eta, sd = sigma[z])
#' dat <- data.frame(y = y, x1 = x1, x2 = x2, x3 = x3)
#'
#' # x1 and x3's effects differ by group (heterogeneous); x2's does not
#' fit <- westMR(y ~ x1 + x2 + x3, data = dat, G_max = 3, procedure = "effects")
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

# ---- steps_table ------------------------------------------------------------

#' Print the step-by-step selection table
#'
#' Generic function that prints, step by step, the candidates considered and
#' chosen at each iteration of a stepwise search. Dispatches on the class of
#' \code{x} (\code{select_variables} or \code{determine_effects}).
#'
#' @param x An object of class \code{determine_effects} or
#'   \code{select_variables}. These are exactly the \code{effect_determination}
#'   and \code{variable_selection} components of the list returned by
#'   \code{\link{westMR}} (present only when \code{procedure} included
#'   \code{"effects"} or \code{"variables"}, respectively).
#' @param ... Currently unused.
#'
#' @return \code{x}, invisibly.
#' @seealso \code{\link{westMR}}, whose \strong{Value} section documents
#'   where \code{variable_selection} and \code{effect_determination} come
#'   from.
#' @export
#'
#' @examples
#' set.seed(1)
#' n <- 500
#' pi <- c(0.2, 0.3, 0.5)
#' z <- sample(1:3, n, replace = TRUE, prob = pi)
#' x1 <- rnorm(n)
#' x2 <- rnorm(n)
#' x3 <- rnorm(n)
#' x4 <- rnorm(n)
#' sigma <- c(1, 0.5, 2)
#' beta <- rbind(
#'   g1 = c(intercept = -4, x1 = -3, x2 = 1.5, x3 = 2, x4 = 0),
#'   g2 = c(intercept = 0, x1 = 1, x2 = 1.5, x3 = -1, x4 = 0),
#'   g3 = c(intercept = 4, x1 = 3, x2 = 1.5, x3 = 0.5, x4 = 0)
#' )
#' eta <- beta[z, "intercept"] + beta[z, "x1"] * x1 + beta[z, "x2"] * x2 +
#'   beta[z, "x3"] * x3 + beta[z, "x4"] * x4
#' y <- rnorm(n, mean = eta, sd = sigma[z])
#' dat <- data.frame(y = y, x1 = x1, x2 = x2, x3 = x3, x4 = x4)
#'
#' fit <- westMR(
#'   y ~ x1 + x2 + x3 + x4, data = dat, G_max = 3, procedure = "variables"
#' )
#' # procedure = "variables" is what populates fit$variable_selection below
#' steps_table(fit$variable_selection)
steps_table <- function(x, ...) UseMethod("steps_table")

#' Print the Step-by-Step Effect-Type Determination Table
#'
#' For each step of the search, prints the current heterogeneous/homogeneous
#' state, every candidate's weighted p-value, and which candidate (if any)
#' was chosen.
#'
#' @param x An object of class \code{determine_effects}. This is produced
#'   internally by \code{determine_effects()}, and surfaces to users as the
#'   \code{effect_determination} component of the list returned by
#'   \code{\link{westMR}} whenever \code{procedure} includes \code{"effects"}
#'   (see \code{westMR}'s \strong{Value} section).
#' @param ... Currently unused.
#'
#' @return \code{x}, invisibly.
#' @seealso \code{\link{westMR}}
#' @export
#'
#' @examples
#' set.seed(1)
#' n <- 500
#' pi <- c(0.2, 0.3, 0.5)
#' z <- sample(1:3, n, replace = TRUE, prob = pi)
#' x1 <- rnorm(n)
#' x2 <- rnorm(n)
#' x3 <- rnorm(n)
#' sigma <- c(1, 0.5, 2)
#' beta <- rbind(
#'   g1 = c(intercept = -4, x1 = -3, x2 = 1.5, x3 = 2),
#'   g2 = c(intercept = 0, x1 = 1, x2 = 1.5, x3 = -1),
#'   g3 = c(intercept = 4, x1 = 3, x2 = 1.5, x3 = 0.5)
#' )
#' eta <- beta[z, "intercept"] + beta[z, "x1"] * x1 + beta[z, "x2"] * x2 +
#'   beta[z, "x3"] * x3
#' y <- rnorm(n, mean = eta, sd = sigma[z])
#' dat <- data.frame(y = y, x1 = x1, x2 = x2, x3 = x3)
#'
#' fit <- westMR(y ~ x1 + x2 + x3, data = dat, G_max = 3, procedure = "effects")
#' # procedure = "effects" is what populates fit$effect_determination below
#' steps_table(fit$effect_determination)
steps_table.determine_effects <- function(x, ...) {
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
  print_steps(
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
#'   \code{\link{westMR}} whenever \code{procedure} includes \code{"variables"}
#'   (see \code{westMR}'s \strong{Value} section).
#' @param ... Currently unused.
#'
#' @return \code{x}, invisibly.
#' @seealso \code{\link{westMR}}
#' @export
#'
#' @examples
#' set.seed(1)
#' n <- 500
#' pi <- c(0.2, 0.3, 0.5)
#' z <- sample(1:3, n, replace = TRUE, prob = pi)
#' x1 <- rnorm(n)
#' x2 <- rnorm(n)
#' x3 <- rnorm(n)
#' x4 <- rnorm(n)
#' sigma <- c(1, 0.5, 2)
#' beta <- rbind(
#'   g1 = c(intercept = -4, x1 = -3, x2 = 1.5, x3 = 2, x4 = 0),
#'   g2 = c(intercept = 0, x1 = 1, x2 = 1.5, x3 = -1, x4 = 0),
#'   g3 = c(intercept = 4, x1 = 3, x2 = 1.5, x3 = 0.5, x4 = 0)
#' )
#' eta <- beta[z, "intercept"] + beta[z, "x1"] * x1 + beta[z, "x2"] * x2 +
#'   beta[z, "x3"] * x3 + beta[z, "x4"] * x4
#' y <- rnorm(n, mean = eta, sd = sigma[z])
#' dat <- data.frame(y = y, x1 = x1, x2 = x2, x3 = x3, x4 = x4)
#'
#' fit <- westMR(
#'   y ~ x1 + x2 + x3 + x4, data = dat, G_max = 3, procedure = "variables"
#' )
#' # procedure = "variables" is what populates fit$variable_selection below
#' steps_table(fit$variable_selection)
steps_table.select_variables <- function(x, ...) {
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
  print_steps(x$steps, x$predictors, x$direction, label_state)
  invisible(x)
}

print_steps <- function(steps, all_predictors, direction, label_state) {
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
#'   \code{\link{westMR}} (always present, regardless of \code{procedure}; see
#'   \code{westMR}'s \strong{Value} section).
#' @param ... Currently unused.
#'
#' @return \code{x}, invisibly.
#' @seealso \code{\link{westMR}}
#' @export
#'
#' @examples
#' set.seed(1)
#' n <- 500
#' pi <- c(0.2, 0.3, 0.5)
#' z <- sample(1:3, n, replace = TRUE, prob = pi)
#' x1 <- rnorm(n)
#' x2 <- rnorm(n)
#' x3 <- rnorm(n)
#' x4 <- rnorm(n)
#' sigma <- c(1, 0.5, 2)
#' beta <- rbind(
#'   g1 = c(intercept = -4, x1 = -3, x2 = 1.5, x3 = 2, x4 = 0),
#'   g2 = c(intercept = 0, x1 = 1, x2 = 1.5, x3 = -1, x4 = 0),
#'   g3 = c(intercept = 4, x1 = 3, x2 = 1.5, x3 = 0.5, x4 = 0)
#' )
#' eta <- beta[z, "intercept"] + beta[z, "x1"] * x1 + beta[z, "x2"] * x2 +
#'   beta[z, "x3"] * x3 + beta[z, "x4"] * x4
#' y <- rnorm(n, mean = eta, sd = sigma[z])
#' dat <- data.frame(y = y, x1 = x1, x2 = x2, x3 = x3, x4 = x4)
#'
#' fit <- westMR(
#'   y ~ x1 + x2 + x3 + x4, data = dat, G_max = 3, procedure = "variables"
#' )
#' # fit$best_fit is always populated, at the BIC-optimal G
#' print(fit$best_fit)
print.fit_fmr <- function(x, ...) {
  cat("Best fit (G: ", x$G, ")\n", sep = "")

  if (!is.null(x$parameter_values$pi_g)) {
    cat("\nMixing proportions:\n")
    print(round(x$parameter_values$pi_g, 4))
  }

  if (!is.null(x$parameter_values$sigma_g)) {
    cat("\nComponent standard deviations:\n")
    print(round(x$parameter_values$sigma_g, 4))
  }

  if (!is.null(x$parameter_values$beta_g)) {
    cat("\nHeterogeneous coefficients:\n")
    print(round(x$parameter_values$beta_g, 4))
  }

  beta_present <- !is.null(x$parameter_values$beta) &&
    length(x$parameter_values$beta) > 0L
  if (beta_present) {
    cat("\nHomogeneous coefficients:\n")
    print(round(x$parameter_values$beta, 4))
  }

  invisible(x)
}

# ---- westMR -----------------------------------------------------------------

#' Print a westMR Result
#'
#' Prints a top-level summary of a \code{\link{westMR}} fit: the family, the
#' candidate \code{G} values, and which procedure(s) were run, followed by the
#' variable-selection result (if run), the effect-determination result (if
#' run), and the best fit.
#'
#' @param x An object of class \code{westMR}, as returned by
#'   \code{\link{westMR}}.
#' @param ... Currently unused.
#'
#' @return \code{x}, invisibly.
#' @export
#'
#' @examples
#' set.seed(1)
#' n <- 500
#' pi <- c(0.2, 0.3, 0.5)
#' z <- sample(1:3, n, replace = TRUE, prob = pi)
#' x1 <- rnorm(n)
#' x2 <- rnorm(n)
#' x3 <- rnorm(n)
#' x4 <- rnorm(n)
#' sigma <- c(1, 0.5, 2)
#' beta <- rbind(
#'   g1 = c(intercept = -4, x1 = -3, x2 = 1.5, x3 = 2, x4 = 0),
#'   g2 = c(intercept = 0, x1 = 1, x2 = 1.5, x3 = -1, x4 = 0),
#'   g3 = c(intercept = 4, x1 = 3, x2 = 1.5, x3 = 0.5, x4 = 0)
#' )
#' eta <- beta[z, "intercept"] + beta[z, "x1"] * x1 + beta[z, "x2"] * x2 +
#'   beta[z, "x3"] * x3 + beta[z, "x4"] * x4
#' y <- rnorm(n, mean = eta, sd = sigma[z])
#' dat <- data.frame(y = y, x1 = x1, x2 = x2, x3 = x3, x4 = x4)
#'
#' fit <- westMR(y ~ x1 + x2 + x3 + x4, data = dat, G_max = 3)
#' print(fit)
print.westMR <- function(x, ...) {
  cat("westMR | family: ", x$family,
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

#' Summarize a westMR Result
#'
#' Collects the call, formula, family, procedure, candidate \code{G} values, best
#' fit, and (if run) the variable-selection and effect-determination results
#' and their step-by-step tables into a single \code{summary.westMR} object
#' for printing.
#'
#' @param object An object of class \code{westMR}, as returned by
#'   \code{\link{westMR}}.
#' @param ... Currently unused.
#'
#' @return An object of class \code{summary.westMR}: a list containing the
#'   call, formula, family, procedure, candidate \code{G} values, the best fit,
#'   the variable-selection and effect-determination results (if run), and
#'   their step-by-step tables.
#' @export
#'
#' @examples
#' set.seed(1)
#' n <- 500
#' pi <- c(0.2, 0.3, 0.5)
#' z <- sample(1:3, n, replace = TRUE, prob = pi)
#' x1 <- rnorm(n)
#' x2 <- rnorm(n)
#' x3 <- rnorm(n)
#' x4 <- rnorm(n)
#' sigma <- c(1, 0.5, 2)
#' beta <- rbind(
#'   g1 = c(intercept = -4, x1 = -3, x2 = 1.5, x3 = 2, x4 = 0),
#'   g2 = c(intercept = 0, x1 = 1, x2 = 1.5, x3 = -1, x4 = 0),
#'   g3 = c(intercept = 4, x1 = 3, x2 = 1.5, x3 = 0.5, x4 = 0)
#' )
#' eta <- beta[z, "intercept"] + beta[z, "x1"] * x1 + beta[z, "x2"] * x2 +
#'   beta[z, "x3"] * x3 + beta[z, "x4"] * x4
#' y <- rnorm(n, mean = eta, sd = sigma[z])
#' dat <- data.frame(y = y, x1 = x1, x2 = x2, x3 = x3, x4 = x4)
#'
#' fit <- westMR(y ~ x1 + x2 + x3 + x4, data = dat, G_max = 3)
#' summary(fit)
summary.westMR <- function(object, ...) {
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
  class(out) <- "summary.westMR"
  out
}

#' Print a westMR Summary
#'
#' Prints the top-level summary produced by \code{\link{summary.westMR}},
#' including the step-by-step tables for the variable-selection and
#' effect-determination searches (if run), and the best fit.
#'
#' @param x An object of class \code{summary.westMR}, as returned by
#'   \code{\link{summary.westMR}}.
#' @param ... Currently unused.
#'
#' @return \code{x}, invisibly.
#' @export
#'
#' @examples
#' set.seed(1)
#' n <- 500
#' pi <- c(0.2, 0.3, 0.5)
#' z <- sample(1:3, n, replace = TRUE, prob = pi)
#' x1 <- rnorm(n)
#' x2 <- rnorm(n)
#' x3 <- rnorm(n)
#' x4 <- rnorm(n)
#' sigma <- c(1, 0.5, 2)
#' beta <- rbind(
#'   g1 = c(intercept = -4, x1 = -3, x2 = 1.5, x3 = 2, x4 = 0),
#'   g2 = c(intercept = 0, x1 = 1, x2 = 1.5, x3 = -1, x4 = 0),
#'   g3 = c(intercept = 4, x1 = 3, x2 = 1.5, x3 = 0.5, x4 = 0)
#' )
#' eta <- beta[z, "intercept"] + beta[z, "x1"] * x1 + beta[z, "x2"] * x2 +
#'   beta[z, "x3"] * x3 + beta[z, "x4"] * x4
#' y <- rnorm(n, mean = eta, sd = sigma[z])
#' dat <- data.frame(y = y, x1 = x1, x2 = x2, x3 = x3, x4 = x4)
#'
#' fit <- westMR(y ~ x1 + x2 + x3 + x4, data = dat, G_max = 3)
#' print(summary(fit))
print.summary.westMR <- function(x, ...) {
  cat("westMR summary | family: ", x$family,
    " | G: ", paste(x$G_values, collapse = ","),
    " | procedure: ", x$procedure, "\n\n",
    sep = ""
  )

  if (!is.null(x$variable_selection)) {
    print(x$variable_selection)
    cat("\n")
    steps_table(x$variable_selection)
  }

  if (!is.null(x$effect_determination)) {
    print(x$effect_determination)
    cat("\n")
    steps_table(x$effect_determination)
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
