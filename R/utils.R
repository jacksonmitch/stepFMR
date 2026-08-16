common_eta <- function(B, beta = numeric(0)) {
  B <- as.matrix(B)
  beta <- as.numeric(beta)

  if (ncol(B) == 0) {
    return(rep(0, nrow(B)))
  }

  if (length(beta) != ncol(B)) {
    stop(
      "Length of beta must equal ncol(B). Got length(beta) = ",
      length(beta),
      " and ncol(B) = ",
      ncol(B),
      "."
    )
  }

  drop(B %*% beta)
}

linear_predictor_matrix <- function(A, B, beta_g, beta) {
  A <- as.matrix(A)

  if (is.null(B)) {
    B <- matrix(numeric(0), nrow = nrow(A), ncol = 0)
  } else {
    B <- as.matrix(B)
  }

  beta_g <- as.matrix(beta_g)
  beta <- as.numeric(beta)

  eta <- (A %*% t(beta_g)) + common_eta(B, beta)
  eta <- as.matrix(eta)
  storage.mode(eta) <- "double"

  eta
}

format_none <- function(x) {
  if (is.null(x) || length(x) == 0L) {
    return("none")
  }

  paste(x, collapse = ", ")
}

coerce_binomial_formula <- function(formula, data) {
  # Just in case
  if (".binom_size" %in% names(data)) {
    stop("Column name '.binom_size' is reserved; please rename your column.")
  }

  response <- formula[[2]]

  if (is.call(response) && identical(response[[1]], as.name("cbind"))) {
    success_expr <- response[[2]]
    failure_expr <- response[[3]]

    successes <- eval(success_expr, data)
    failures <- eval(failure_expr, data)

    if (any(successes < 0) || any(failures < 0)) {
      stop("Successes and failures must be non-negative.")
    }

    data$.binom_size <- successes + failures
    formula[[2]] <- success_expr # response becomes the success count only instead of am ugly cbind

    return(list(formula = formula, data = data))
  }

  # Plain response: treat as 0/1 (Bernoulli, size = 1)
  y <- eval(response, data)

  if (!all(y %in% c(0, 1))) {
    stop("Binomial response must be cbind(successes, failures) or a 0/1 vector.")
  }

  data$.binom_size <- rep(1, length(y))
  warning(
    "Received Binomial family with Bernoulli/binary response. ",
    "Each observation carries very little information about its group ",
    "membership, so mixture regression fits with binary responses tend to ",
    "be unstable and imprecise. Grouped binomial counts (e.g. cbind(successes, ",
    "failures) ~ .) are recommended when available.",
    call. = FALSE
  )
  list(formula = formula, data = data)
}

compact <- function(x) Filter(Negate(is.null), x)
