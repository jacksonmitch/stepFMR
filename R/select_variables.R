# select_variables.R

#' Stepwise Selection of Predictors
#'
#' Runs the WEST stepwise procedure to decide which predictors should be
#' included in the model. At each step, remaining candidates are tested
#' against a shared baseline fit via \code{test_predictors()}, and the most
#' eligible candidate (if any) is added (forward) or removed (backward).
#'
#' @param model A \code{FMRModel} object.
#' @param direction A character string, either \code{"forward"} or
#'   \code{"backward"}, specifying the stepwise search direction.
#' @param predictors A character vector of candidate predictor names to
#'   consider for selection. Defaults to all predictors in \code{model}.
#'
#' @return A list of class \code{select_variables} with elements:
#'   \code{direction}, \code{alpha}, \code{predictors} (all candidate
#'   predictor names considered), \code{selected} (the chosen predictor
#'   names), \code{steps} (the step-by-step search log), \code{final_fits}
#'   (fits across \code{model$G_values} for the final selected model), and
#'   \code{final_formula}.
#' @noRd
select_variables <- function(
  model,
  direction = "forward",
  predictors = model$predictors
) {
  stopifnot(direction %in% c("forward", "backward"))

  if (direction == "forward") {
    included <- character(0)
    update_included <- union
  } else {
    included <- predictors
    update_included <- setdiff
  }

  result <- test_predictors(
    model = model,
    predictors = predictors,
    included = included,
    heterogeneous = included, # selection: het always mirrors included
    update_included = update_included,
    direction = direction
  )

  out <- list(
    direction = direction,
    alpha = model$control$alpha,
    predictors = predictors,
    selected = result$included,
    steps = result$steps,
    final_fits = result$final_fits,
    final_formula = make_formula(result$included, response = model$response)
  )

  class(out) <- "select_variables"
  out
}
