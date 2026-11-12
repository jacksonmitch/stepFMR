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
#   "kmeans_2" -> "kmeans", "random_4" -> "random"
#   "quantile_predictor_x1" -> "quantile_predictor"
#   "previous_fit" (warm start) and "quantile_response" are unchanged.
#
# Only the n_best_init starts that survive burn-in are run to convergence, so
# "final" comparisons only cover those. To compare every start after full EM,
# run_init_study() defaults to keeping ALL starts (slower), including extra
# deterministic and warm starts beyond n_init, to evaluate smaller cutoffs.
# These are counterfactual comparisons on the recorded search path; changing
# n_best_init can change the path of a fresh stepFMR() run.

rm(list = ls())
devtools::load_all()

# ---- helpers ----------------------------------------------------------------

init_type <- function(name) {
  known <- c(
    "kmeans_raw", "kmeans", "mclust", "random", "previous_fit",
    "quantile_residual", "quantile_response", "quantile_predictor",
    "quantile_neg_A", "quantile_A"
  )
  vapply(name, function(nm) {
    hit <- known[startsWith(nm, known)]
    if (length(hit) == 0L) sub("_*[0-9]+$", "", nm) else hit[1]
  }, character(1), USE.NAMES = FALSE)
}

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
extract_fit <- function(fit, stage, step, candidate, tie_tol = 0.01) {
  init <- fit[["initialization"]]
  final_ll <- finite_sorted(init[["final_logliks"]])
  burnin_ll <- finite_sorted(init[["initial_logliks"]])
  winner <- fit[["best_init_name"]]
  if (length(final_ll) == 0L || length(burnin_ll) == 0L) return(NULL)

  # sort() preserves input order for ties, as does select_best_initialization().
  # Use positions rather than rank(): tied starts still occupy separate slots.
  winner_rank <- match(winner, names(burnin_ll))
  final_by_rank <- unname(init[["final_logliks"]][names(burnin_ll)])
  competitive <- which(is.finite(final_by_rank) &
    final_by_rank >= max(final_ll) - tie_tol)

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
    winner_burnin_rank = winner_rank,
    competitive_burnin_rank = if (length(competitive)) min(competitive) else NA_integer_,
    init_burnin = if (is.null(init[["burnin"]])) NA_integer_ else init[["burnin"]],
    n_burnin_valid = length(burnin_ll),
    n_final = sum(is.finite(final_by_rank)),
    all_starts_evaluated = all(is.finite(final_by_rank)),
    stringsAsFactors = FALSE
  ))

  type_rows <- rbind(
    cbind(where, basis = "final", final_types),
    cbind(where, basis = "burnin", strategy_table(burnin_ll))
  )

  best_by_k <- cummax(ifelse(is.finite(final_by_rank), final_by_rank, -Inf))
  top_k <- cbind(where, data.frame(
    n_best_init = seq_along(burnin_ll),
    n_burnin_valid = length(burnin_ll),
    winner_kept = seq_along(burnin_ll) >= winner_rank,
    loss = ifelse(is.finite(best_by_k), pmax(0, max(final_ll) - best_by_k), NA_real_)
  ))

  list(fit = fit_row, types = type_rows, top_k = top_k)
}

# Walk every fit stored in one stage's step log (variable_selection or
# effect_determination)
collect_stage <- function(stage_result, stage, tie_tol = 0.01) {
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
        fit, stage, job[["step"]], job[["candidate"]], tie_tol
      )
    }
  }

  list(
    fits = do.call(rbind, lapply(extracted, `[[`, "fit")),
    types = do.call(rbind, lapply(extracted, `[[`, "types")),
    top_k = do.call(rbind, lapply(extracted, `[[`, "top_k"))
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

# Burn-in position of the exact winner (not its initialization name's suffix).
burnin_rank_table <- function(fits, by = c("stage", "context")) {
  groups <- split(fits, fits[by], drop = TRUE)
  out <- do.call(rbind, lapply(groups, function(d) {
    counts <- table(d[["winner_burnin_rank"]], useNA = "ifany")
    cbind(d[rep(1L, length(counts)), by, drop = FALSE], data.frame(
      winner_burnin_rank = as.integer(names(counts)),
      fits = as.integer(counts),
      pct = 100 * as.integer(counts) / nrow(d)
    ))
  }))
  rownames(out) <- NULL
  out
}

# One row per cutoff, pooled over fits. For fits with fewer valid starts than k,
# keeping k means keeping all available starts. Missing final logliks cannot
# establish the quality of an untested start; losses use observed finals only.
retention_table <- function(top_k, tie_tol, by = c("stage", "context")) {
  groups <- split(top_k, top_k[by], drop = TRUE)
  out <- do.call(rbind, lapply(groups, function(group) {
    do.call(rbind, lapply(seq_len(max(group[["n_best_init"]])), function(k) {
      d <- group[group[["n_best_init"]] == pmin(k, group[["n_burnin_valid"]]), ]
      known <- !is.na(d[["loss"]])
      missed <- known & d[["loss"]] > tie_tol
      row <- cbind(d[1, by, drop = FALSE], data.frame(
        n_best_init = k,
        fits = nrow(d),
        winner_dropped = sum(!d[["winner_kept"]], na.rm = TRUE),
        pct_winner_kept = 100 * mean(d[["winner_kept"]], na.rm = TRUE),
        unknown_loss = sum(!known),
        loss_gt_tol = sum(missed),
        pct_within_tol = if (any(known)) 100 * mean(d[["loss"]][known] <= tie_tol) else NA_real_,
        mean_loss = or_na(d[["loss"]], mean),
        max_loss = or_na(d[["loss"]], max)
      ))
      if ("seed" %in% names(d)) row[["seeds_with_loss"]] <- length(unique(d[["seed"]][missed]))
      row
    }))
  }))
  rownames(out) <- NULL
  out
}

print_burnin_diagnostics <- function(x) {
  fits <- x[["fits"]]
  if (is.null(fits[["winner_burnin_rank"]])) {
    cat("\nBurn-in ranks unavailable in this saved summary; rerun the study.\n")
    return(invisible(NULL))
  }
  cat("\n== Final winner's rank at the end of burn-in ==\n")
  cat("Rank 1 is the burn-in leader; exact ties retain their selection order.\n")
  print(x[["burnin_ranks"]], row.names = FALSE, digits = 3)
  cat(sprintf("\nAll valid starts evaluated to full EM: %d / %d fits\n",
    sum(fits[["all_starts_evaluated"]]), nrow(fits)))
  if (!all(fits[["all_starts_evaluated"]])) {
    cat("Partial coverage: ranks/losses compare observed finals only; discarded starts may win.\n")
  }
  cat("Burn-in iterations:", paste(sort(unique(fits[["init_burnin"]])), collapse = ", "), "\n")
  cat(sprintf("Smallest k retaining every observed winner: %g; within tie_tol: %g\n",
    or_na(fits[["winner_burnin_rank"]], max),
    or_na(fits[["competitive_burnin_rank"]], max)))
  cat("\n== What if only the top k burn-in starts were kept? ==\n")
  cat("winner_dropped: exact winner excluded; loss_gt_tol: best retained final loses > tie_tol.\n")
  cat("pct_within_tol excludes unknown_loss fits (no finite final among retained starts).\n")
  cat("These comparisons use the recorded search path; rerun a chosen cutoff to confirm.\n")
  print(x[["retention"]], row.names = FALSE, digits = 3)
  cat("\nFits needing >2 starts for a final within tie_tol: ",
    sum(fits[["competitive_burnin_rank"]] > 2L, na.rm = TRUE),
    " (inspect $fits, or $top_k for per-fit losses)\n", sep = "")
  invisible(x)
}

# ---- main -------------------------------------------------------------------

# `result` is the object returned by stepFMR(); running this on a saved result
# avoids refitting when only the summaries change.
init_diagnostics <- function(result, tie_tol = 0.01) {
  stopifnot(inherits(result, "stepFMR"))
  stopifnot(length(tie_tol) == 1L, is.finite(tie_tol), tie_tol >= 0)

  stages <- compact(list(
    variables = result[["variable_selection"]],
    effects = result[["effect_determination"]]
  ))
  collected <- lapply(names(stages), function(nm) {
    collect_stage(stages[[nm]], nm, tie_tol)
  })

  fits <- do.call(rbind, lapply(collected, `[[`, "fits"))
  types <- do.call(rbind, lapply(collected, `[[`, "types"))
  top_k <- do.call(rbind, lapply(collected, `[[`, "top_k"))
  if (is.null(fits)) {
    stop("No search steps found in `result`; nothing to summarise.")
  }

  fits[["context"]] <- ifelse(fits[["step"]] == 0L, "baseline", "candidate")
  types[["context"]] <- ifelse(types[["step"]] == 0L, "baseline", "candidate")
  top_k[["context"]] <- ifelse(top_k[["step"]] == 0L, "baseline", "candidate")

  out <- list(
    tie_tol = tie_tol,
    fits = fits,
    types = types,
    top_k = top_k,
    burnin_ranks = burnin_rank_table(fits),
    retention = retention_table(top_k, tie_tol),
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
  print_burnin_diagnostics(x)

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

# ---- multi-seed study -------------------------------------------------------

# Aggregate a named list of init_diagnostics() results (one per seed)
aggregate_init_study <- function(studies, tie_tol = 0.01) {
  seeds <- names(studies)

  types <- do.call(rbind, lapply(seeds, function(s) {
    cbind(seed = s, studies[[s]][["types"]])
  }))
  fits <- do.call(rbind, lapply(seeds, function(s) {
    cbind(seed = s, studies[[s]][["fits"]])
  }))
  top_k <- do.call(rbind, lapply(seeds, function(s) {
    if (is.null(studies[[s]][["top_k"]])) {
      stop("Study for seed ", s, " lacks burn-in ranks; rerun init_diagnostics() on its stepFMR result.")
    }
    cbind(seed = s, studies[[s]][["top_k"]])
  }))

  types[["needed"]] <- !is.na(types[["margin"]]) & types[["margin"]] > tie_tol

  # Per seed, then pooled across seeds
  per_seed <- do.call(rbind, lapply(
    split(types, types[c("seed", "stage", "context", "basis", "type")], drop = TRUE),
    function(d) {
      data.frame(
        seed = d[["seed"]][1],
        stage = d[["stage"]][1],
        context = d[["context"]][1],
        basis = d[["basis"]][1],
        type = d[["type"]][1],
        fits = nrow(d),
        needed = sum(d[["needed"]]),
        loss_sum = sum(d[["margin"]][d[["needed"]]]),
        loss_max = or_na(d[["margin"]], max),
        stringsAsFactors = FALSE
      )
    }
  ))

  marginal <- do.call(rbind, lapply(
    split(per_seed, per_seed[c("stage", "context", "basis", "type")], drop = TRUE),
    function(d) {
      pct_by_seed <- 100 * d[["needed"]] / d[["fits"]]
      data.frame(
        stage = d[["stage"]][1],
        context = d[["context"]][1],
        basis = d[["basis"]][1],
        type = d[["type"]][1],
        fits = sum(d[["fits"]]),
        of = sum(fits[["stage"]] == d[["stage"]][1] &
          fits[["context"]] == d[["context"]][1]),
        needed = sum(d[["needed"]]),
        pct_needed = 100 * sum(d[["needed"]]) / sum(d[["fits"]]),
        pct_seed_min = min(pct_by_seed),
        pct_seed_max = max(pct_by_seed),
        seeds_needed = sum(d[["needed"]] > 0),
        mean_loss_when_needed = if (sum(d[["needed"]]) > 0) {
          sum(d[["loss_sum"]]) / sum(d[["needed"]])
        } else {
          NA_real_
        },
        max_loss = or_na(d[["loss_max"]], max),
        stringsAsFactors = FALSE
      )
    }
  ))
  marginal <- marginal[
    order(
      marginal[["stage"]], marginal[["context"]],
      marginal[["basis"]], -marginal[["needed"]]
    ),
  ]
  rownames(marginal) <- NULL

  out <- list(
    seeds = seeds,
    tie_tol = tie_tol,
    fits = fits,
    types = types,
    top_k = top_k,
    burnin_ranks = burnin_rank_table(fits),
    burnin_ranks_by_seed = burnin_rank_table(fits, c("seed", "stage", "context")),
    retention = retention_table(top_k, tie_tol),
    retention_by_seed = retention_table(top_k, tie_tol, c("seed", "stage", "context")),
    wins_by_step = win_table(fits, c("stage", "step"), tie_tol),
    wins_overall = win_table(fits, "context", tie_tol),
    marginal = marginal
  )
  class(out) <- "init_study"
  out
}

print.init_study <- function(x, ...) {
  cat(sprintf(
    "Seeds: %d  |  fits: %d  |  tie tolerance: %g  |  burn-in leader == final winner: %.1f%%\n",
    length(x[["seeds"]]), nrow(x[["fits"]]), x[["tie_tol"]],
    100 * mean(x[["fits"]][["burnin_agrees"]])
  ))
  print_burnin_diagnostics(x)

  cat("\n== Loss in best loglik if a strategy were dropped ==\n")
  cat("fits: fits where the strategy ran (finite), of: all fits in the context\n")
  cat("pct_needed: needed / fits; pct_seed_min/max: range across seeds\n")
  cat("seeds_needed: seeds where it was needed at least once\n")
  print(x[["marginal"]], row.names = FALSE, digits = 3)

  cat("\n== Winners overall ==\n")
  print(x[["wins_overall"]], row.names = FALSE, digits = 3)

  cat("\n(per-step winners: $wins_by_step; per-fit detail: $fits, $types)\n")
  invisible(x)
}

run_init_study <- function(seeds,
                           params,
                           n = 500,
                           G_max = 4,
                           control = fmr_control(parallel = TRUE),
                           tie_tol = 0.1,
                           keep_all = TRUE) {
  # Diagnostic-only override: n_init is a minimum pool size, so deterministic
  # starts and previous_fit can push the actual count above it. The selector
  # caps n_keep at the number of valid starts. Leave pool generation unchanged.
  if (keep_all) control[["n_best_init"]] <- .Machine$integer.max
  if (keep_all) cat("Reference study: running full EM from every valid burn-in start.\n")
  # Note: binomial with size > 1 needs a cbind(success, failure) response,
  # which simulate_fmr() does not produce
  family <- if (is.null(params[["family"]])) "gaussian" else params[["family"]]

  studies <- list()
  errors <- list()

  for (seed in seeds) {
    cat(sprintf("Seed %d (%d of %d)\n", seed, match(seed, seeds), length(seeds)))

    one <- tryCatch(
      {
        sim <- do.call(simulate_fmr, c(params, list(n = n, seed = seed)))
        result <- stepFMR(
          sim$formula, sim$data,
          G_max = G_max,
          family = family,
          control = control
        )
        init_diagnostics(result, tie_tol = tie_tol)
      },
      error = function(e) e
    )

    if (inherits(one, "error")) {
      errors[[as.character(seed)]] <- conditionMessage(one)
    } else {
      studies[[as.character(seed)]] <- one
    }
  }

  if (length(studies) == 0L) {
    stop("Every seed failed. First error: ", errors[[1]])
  }
  if (length(errors) > 0L) {
    warning(sprintf("%d of %d seeds failed; see $errors", length(errors), length(seeds)))
  }

  out <- aggregate_init_study(studies, tie_tol)
  out[["errors"]] <- errors
  out[["keep_all"]] <- keep_all
  out
}
# ---- run --------------------------------------------------------------------

# Swap in any dataset/formula here (binomial data needs a .binom_size column)

sim <- do.call(
  simulate_fmr,
  c(scenarios$four_group_twelve_variables_poisson, list(n = 300, seed = 2))
)

# Remember to devtools::install() so parallel workers see current code
result <- stepFMR(
  sim$formula, sim$data,
  G_max = 3,
  family = "poisson",
  direction = "backward",
  control = fmr_control(parallel = FALSE, n_best_init = 2, init_burnin = 30, irwls_max_iter = 2 )
)

init_diag_poisson_3 <- init_diagnostics(result, tie_tol = 1)
print(init_diag_poisson_3)

# Multi-seed usage: keep_all = TRUE (the default) supplies a full reference for
# choosing n_best_init; use FALSE to diagnose only the configured survivors.
# Compare runs with different init_burnin values to tune burn-in length too.
# Useful details after each run:
#   study$retention                       # loss for each cutoff, including k=2
#   study$retention_by_seed               # detect seed-specific failures
#   subset(study$fits, winner_burnin_rank > 2)
#   subset(study$fits, competitive_burnin_rank > 2) # material top-2 misses
NoSleepR::with_nosleep({
  study <- run_init_study(1:100, 
    scenarios$four_group_twelve_variables_binomial,
    n = 500,
    G_max = 5,
    control = fmr_control(n_best_init = 10, init_burnin = 50, verbose = FALSE),
    tie_tol = 1
  )
})
print(study)
saveRDS(study, "init_study_binomial_50burnin_allstarts.rds")

NoSleepR::with_nosleep({
  study <- run_init_study(1:100, 
    scenarios$four_group_twelve_variables_poisson,
    n = 500,
    G_max = 5,
    control = fmr_control(n_best_init = 10, init_burnin = 50, verbose = FALSE),
    tie_tol = 1
  )
})
print(study)
saveRDS(study, "init_study_poisson_50burnin_allstarts.rds")

NoSleepR::with_nosleep({
  study <- run_init_study(1:100, 
    scenarios$four_group_twelve_variables,
    n = 500,
    G_max = 5,
    control = fmr_control(n_best_init = 10, init_burnin = 50, verbose = FALSE),
    tie_tol = 1
  )
})
print(study)
saveRDS(study, "init_study_gaussian_50burnin_allstarts.rds")
