# Initialization diagnostics
#
# Runs stepFMR on a dataset, then tallies which initialization strategy won at
# every fit made during the search (the cold-start baseline fits, plus every
# candidate fit at every step), and by how much.
#
#   margin   winner minus the best start from any OTHER strategy, i.e. what the
#            fit loses if the winner's strategy is removed from the pool.
#
# Starts are grouped into strategies by dropping the numeric suffix:
#   "kmeans_2" -> "kmeans", "random_4" -> "random", "quantile_A3" -> "quantile_A"
#   "previous_fit" (warm start) and "quantile_response" are unchanged.
#
# Only the n_best_init starts that survive burn-in are run to convergence, so
# "final" comparisons only cover those. To compare every start after full EM,
# use build_control(n_best_init = n_init) (slower).

rm(list = ls())
devtools::load_all()

# ---- helpers ----------------------------------------------------------------

init_type <- function(name) sub("_?[0-9]+$", "", name)

finite_sorted <- function(ll) sort(ll[is.finite(ll)], decreasing = TRUE)

# Apply f to x after dropping NAs; NA if nothing is left
or_na <- function(x, f) {
  x <- x[!is.na(x)]
  if (length(x) == 0L) NA_real_ else f(x)
}

# One row per strategy: loglik of its best start, and how far the overall best
# loglik falls if the strategy is removed (NA if it is the only strategy)
strategy_table <- function(ll) {
  types <- init_type(names(ll))
  unique_types <- unique(types)

  data.frame(
    type = unique_types,
    best_ll = vapply(
      unique_types,
      function(t) max(ll[types == t]),
      numeric(1),
      USE.NAMES = FALSE
    ),
    margin = vapply(
      unique_types,
      function(t) {
        others <- ll[types != t]
        if (length(others) == 0L) NA_real_ else max(ll) - max(others)
      },
      numeric(1),
      USE.NAMES = FALSE
    ),
    stringsAsFactors = FALSE
  )
}

# One fit_fmr -> one row in `fit`, and one row per strategy in `types`
extract_fit <- function(fit, stage, step, candidate) {
  init <- fit[["initialization"]]
  final_ll <- finite_sorted(init[["final_logliks"]])
  burnin_ll <- finite_sorted(init[["initial_logliks"]])
  winner <- fit[["best_init_name"]]

  where <- data.frame(
    stage = stage,
    step = step,
    candidate = candidate,
    G = fit[["G"]],
    stringsAsFactors = FALSE
  )

  final_types <- strategy_table(final_ll)
  has_runner_up <- length(final_ll) >= 2L

  fit_row <- cbind(where, data.frame(
    winner = winner,
    winner_type = init_type(winner),
    runner_up = if (has_runner_up) names(final_ll)[2] else NA_character_,
    margin = final_types[["margin"]][final_types[["type"]] == init_type(winner)],
    burnin_top = names(burnin_ll)[1],
    burnin_agrees = identical(names(burnin_ll)[1], winner),
    stringsAsFactors = FALSE
  ))

  type_rows <- rbind(
    cbind(where, basis = "final", final_types),
    cbind(where, basis = "burnin", strategy_table(burnin_ll))
  )

  list(fit = fit_row, types = type_rows)
}

# Walk every fit stored in one stage's step log (variable_selection or
# effect_determination)
collect_stage <- function(stage_result, stage) {
  steps <- stage_result[["steps"]]
  if (length(steps) == 0L) {
    return(NULL)
  }

  # Step 0: the cold-start baseline that the step-1 tests were compared against
  jobs <- list(list(
    fits = steps[[1]][["tests"]][[1]][["shared_fits"]],
    step = 0L,
    candidate = "(baseline)"
  ))
  for (s in steps) {
    for (cand in names(s[["tests"]])) {
      jobs[[length(jobs) + 1L]] <- list(
        fits = s[["tests"]][[cand]][["candidate_fits"]],
        step = s[["step"]],
        candidate = cand
      )
    }
  }

  extracted <- list()
  for (job in jobs) {
    for (fit in job[["fits"]]) {
      extracted[[length(extracted) + 1L]] <- extract_fit(
        fit, stage, job[["step"]], job[["candidate"]]
      )
    }
  }

  list(
    fits = do.call(rbind, lapply(extracted, `[[`, "fit")),
    types = do.call(rbind, lapply(extracted, `[[`, "types"))
  )
}

# Winner counts grouped by `by` (e.g. c("stage", "step")), with how decisive
# the wins were. A win is "decisive" when margin > tie_tol.
win_table <- function(fits, by, tie_tol) {
  fits[["decisive"]] <- !is.na(fits[["margin"]]) & fits[["margin"]] > tie_tol

  groups <- split(fits, fits[c(by, "winner_type")], drop = TRUE)
  rows <- lapply(groups, function(d) {
    cbind(
      d[1, c(by, "winner_type"), drop = FALSE],
      data.frame(
        wins = nrow(d),
        decisive = sum(d[["decisive"]]),
        mean_margin = or_na(d[["margin"]][d[["decisive"]]], mean),
        max_margin = or_na(d[["margin"]], max)
      )
    )
  })

  out <- do.call(rbind, rows)
  out[["share"]] <- out[["wins"]] / stats::ave(out[["wins"]], out[by], FUN = sum)

  out <- out[do.call(order, c(unname(as.list(out[by])), list(-out[["wins"]]))), ]
  out <- out[c(by, "winner_type", "wins", "share", "decisive",
               "mean_margin", "max_margin")]
  rownames(out) <- NULL
  out
}

# For each strategy: in how many fits would dropping it have cost loglik, and
# how much
marginal_table <- function(types, tie_tol) {
  groups <- split(types, types[c("context", "basis", "type")], drop = TRUE)
  rows <- lapply(groups, function(d) {
    data.frame(
      context = d[["context"]][1],
      basis = d[["basis"]][1],
      type = d[["type"]][1],
      fits = nrow(d),
      needed = sum(d[["margin"]] > tie_tol, na.rm = TRUE),
      mean_loss_if_dropped = or_na(d[["margin"]], mean),
      max_loss_if_dropped = or_na(d[["margin"]], max),
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, rows)
  out[["pct_needed"]] <- 100 * out[["needed"]] / out[["fits"]]
  out <- out[order(out[["context"]], out[["basis"]], -out[["needed"]]), ]
  rownames(out) <- NULL
  out
}

# ---- main -------------------------------------------------------------------

# `result` is the object returned by stepFMR(); running this on a saved result
# avoids refitting when only the summaries change.
init_diagnostics <- function(result, tie_tol = 0.01) {
  stopifnot(inherits(result, "stepFMR"))

  stages <- compact(list(
    variables = result[["variable_selection"]],
    effects = result[["effect_determination"]]
  ))
  collected <- lapply(names(stages), function(nm) {
    collect_stage(stages[[nm]], nm)
  })

  fits <- do.call(rbind, lapply(collected, `[[`, "fits"))
  types <- do.call(rbind, lapply(collected, `[[`, "types"))
  if (is.null(fits)) {
    stop("No search steps found in `result`; nothing to summarise.")
  }

  fits[["context"]] <- ifelse(fits[["step"]] == 0L, "baseline", "candidate")
  types[["context"]] <- ifelse(types[["step"]] == 0L, "baseline", "candidate")

  out <- list(
    tie_tol = tie_tol,
    fits = fits,
    types = types,
    wins_by_step = win_table(fits, c("stage", "step"), tie_tol),
    wins_overall = win_table(fits, "context", tie_tol),
    marginal = marginal_table(types, tie_tol)
  )
  class(out) <- "init_diagnostics"
  out
}

print.init_diagnostics <- function(x, ...) {
  fits <- x[["fits"]]
  cat(sprintf(
    "Fits examined: %d  |  tie tolerance: %g loglik units\n",
    nrow(fits), x[["tie_tol"]]
  ))
  cat(sprintf(
    "Burn-in leader was the final winner in %.1f%% of fits\n",
    100 * mean(fits[["burnin_agrees"]])
  ))

  cat("\n== Winners by step ==\n")
  by_step <- x[["wins_by_step"]]
  keys <- unique(by_step[c("stage", "step")])
  for (i in seq_len(nrow(keys))) {
    in_step <- by_step[["stage"]] == keys[["stage"]][i] &
      by_step[["step"]] == keys[["step"]][i]
    rows <- by_step[in_step, setdiff(names(by_step), c("stage", "step"))]

    label <- if (keys[["step"]][i] == 0L) {
      "baseline"
    } else {
      paste("step", keys[["step"]][i])
    }
    cat(sprintf("\n[%s, %s] %d fits\n", keys[["stage"]][i], label, sum(rows[["wins"]])))
    print(rows, row.names = FALSE, digits = 3)
  }

  cat("\n== Winners overall ==\n")
  print(x[["wins_overall"]], row.names = FALSE, digits = 3)

  cat("\n== Loss in best loglik if a strategy were dropped ==\n")
  cat("basis 'final': full-EM logliks of the starts kept after burn-in\n")
  cat("basis 'burnin': burn-in logliks of every start\n")
  print(x[["marginal"]], row.names = FALSE, digits = 3)

  invisible(x)
}

# ---- run --------------------------------------------------------------------

# Swap in any dataset/formula here (binomial data needs a .binom_size column)
sim <- do.call(
  simulate_fmr,
  c(scenarios$four_group_twelve_variables_doubled, list(n = 500, seed = 123))
)

# Remember to devtools::install() so parallel workers see current code
result <- stepFMR(
  sim$formula, sim$data,
  G_max = 4,
  control = fmr_control(parallel = TRUE, n_best_init = 10)
)

init_diag <- init_diagnostics(result, tie_tol = 0.01)
print(init_diag)

# Per-fit detail for digging in:
# View(init_diag$fits)