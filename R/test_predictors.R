# test_predictors.R

#' Shared Stepwise Search Driver
#'
#' Drives the forward/backward stepwise search shared by
#' \code{select_variables()} and \code{determine_effects()}. At each step,
#' it fits a baseline ("shared") model plus one candidate model per
#' remaining covariate, compares each candidate to the shared fit via
#' \code{west_procedure()}, and either adopts the most eligible candidate or
#' stops.
#'
#' @param model A \code{FMRModel} object.
#' @param predictors A character vector of all candidate covariate names
#'   under consideration.
#' @param included A character vector of covariate names currently included
#'   in the model.
#' @param heterogeneous A character vector of covariate names currently
#'   treated as heterogeneous (component-specific).
#' @param update_included A function of \code{(set, predictor)} used to
#'   update the included-covariate set when a candidate is chosen (e.g.
#'   \code{union}/\code{setdiff} for variable selection, or a no-op for
#'   effect-type determination, where the included set never changes).
#' @param direction A character string, either \code{"forward"} or
#'   \code{"backward"}, specifying the stepwise search direction.
#'
#' @return A list with elements: \code{included} (final included covariate
#'   names), \code{heterogeneous} (final heterogeneous covariate names),
#'   \code{common} (\code{included} minus \code{heterogeneous}),
#'   \code{steps} (the step-by-step search log), and \code{final_fits}
#'   (fits across \code{G_values} for the final model).
#' @noRd
test_predictors <- function(
  model,
  predictors,
  included,
  heterogeneous,
  update_included,
  direction = c("forward", "backward"),
  label = "Testing predictors"
) {
  direction <- match.arg(direction)
  alpha <- model$control$alpha

  included <- as.character(included)
  heterogeneous <- as.character(heterogeneous)
  predictors <- as.character(predictors)

  heterogeneous <- intersect(heterogeneous, included)

  # These three always co-vary with direction, not task-specific
  if (direction == "forward") {
    is_eligible <- function(p) is.finite(p) & p < alpha
    find_best <- which.min
    update_heterogeneous <- union
    remaining_predictors <- setdiff(predictors, heterogeneous)
  } else {
    is_eligible <- function(p) is.finite(p) & p >= alpha
    find_best <- which.max
    update_heterogeneous <- setdiff
    remaining_predictors <- intersect(predictors, heterogeneous)
  }

  steps <- list()
  step_id <- 1L

  common <- setdiff(included, heterogeneous)
  shared_data <- prepare_data(model, included = included, common = common)

  progress <- isTRUE(model$control$verbose) && length(remaining_predictors) > 0L
  if (progress) {
    pb <- cli::cli_progress_bar(
      name = label,
      total = length(remaining_predictors),
      current = FALSE,
      auto_terminate = FALSE,
      .auto_close = FALSE,
      format = paste("{trimws(cli::pb_name)}",
                     "{cli::pb_current}/{cli::pb_total} {cli::pb_status}"),
      status = "Fitting initial model"
    )
    on.exit(cli::cli_progress_done(id = pb, result = "auto"), add = TRUE)
    cli::cli_progress_update(id = pb, set = 0, force = TRUE)
  }

  n_best_init <- model$control$n_best_init
  # Always run baseline fit with all inits
  model$control$n_best_init <- model$control$n_init
  shared_fits <- fit_across_G(model, shared_data)
  model$control$n_best_init <- n_best_init
  while (length(remaining_predictors) > 0) {
    tests <- future.apply::future_lapply(remaining_predictors,
      function(predictor) {
        cand_included <- update_included(included, predictor)
        cand_heterogeneous <- update_heterogeneous(heterogeneous, predictor)
        cand_common <- setdiff(cand_included, cand_heterogeneous)

        candidate_data <- prepare_data(model,
          included = cand_included,
          common = cand_common
        )

        west_procedure(
          model = model,
          direction = direction,
          shared_fits = shared_fits,
          candidate_data = candidate_data
        )
      }, future.seed = TRUE)
    names(tests) <- remaining_predictors

    p0 <- vapply(tests, function(t) t$p0, numeric(1))

    steps[[step_id]] <- list(
      step = step_id,
      p0 = p0,
      chosen = NA_character_,
      tests = tests
    )

    eligible_p <- p0[is_eligible(p0)]
    if (length(eligible_p) > 0) {
      chosen <- names(eligible_p)[find_best(eligible_p)]
      steps[[step_id]]$chosen <- chosen
    } else {
      if (progress) {
        cli::cli_progress_update(id = pb, set = step_id,
                                status = "no eligible predictor")
      }
      break
    }

    shared_fits <- tests[[chosen]]$candidate_fits
    included <- update_included(included, chosen)
    heterogeneous <- update_heterogeneous(heterogeneous, chosen)
    remaining_predictors <- setdiff(remaining_predictors, chosen)
    if (progress) {
      cli::cli_progress_update(id = pb, set = step_id,
        status = if (label == "Determining effects:") {
          sprintf("Heterogeneous: %s | Homogeneous: %s",
            cli::ansi_strtrim(format_none(heterogeneous), 40, ellipsis = "..."),
            cli::ansi_strtrim(format_none(setdiff(included, heterogeneous)),
                              30, ellipsis = "...")
          )
        } else if (direction == "backward") {
          sprintf("Excluded: %s", cli::ansi_strtrim(
            format_none(setdiff(predictors, included)), 40, ellipsis = "..."
          ))
        } else {
          sprintf("Included: %s", cli::ansi_strtrim(
            format_none(included), 40, ellipsis = "..."
          ))
        }
      )
    }
    step_id <- step_id + 1L
  }

  list(
    included = included,
    heterogeneous = heterogeneous,
    common = setdiff(included, heterogeneous),
    steps = steps,
    final_fits = shared_fits
  )
}
