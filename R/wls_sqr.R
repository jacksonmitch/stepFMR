#' Weighted Least Squares via Structured QR
#'
#' Solves the stacked weighted least squares problem for a mixture
#' regression with component-specific (\code{X_het}) and shared
#' (\code{X_com}) coefficients without forming normal equations. Stage 1
#' takes a QR decomposition of the weighted heterogeneous design within
#' each component, which eliminates \code{beta_g} from the problem. Stage 2
#' takes a QR decomposition of the stacked residualized common design to
#' solve for the shared \code{beta}. The component coefficients are then
#' recovered by back-substitution.
#'
#' @param dat A \code{FMRData} object.
#' @param z A numeric vector (recycled across components, e.g. the Gaussian
#'   response) or an n x G matrix of working responses.
#' @param weights An n x G matrix of non-negative observation weights.
#' @param control A \code{FMRControl} object. If
#'   \code{return_qr_parts} is \code{TRUE}, the intermediate factors are
#'   also returned.
#'
#' @return A list with \code{beta_g} (G x p_het matrix) and \code{beta}
#'   (numeric vector of length p_com), plus \code{R1_list}, \code{E_list}
#'   and \code{R2} when \code{control$return_qr_parts} is \code{TRUE}.
#' @noRd
wls_sqr <- function(dat, z, weights, control) {
  A <- dat$X_het
  B <- dat$X_com

  stopifnot(!is.null(B))

  w <- weights

  n <- dat$n
  G <- ncol(w)
  p1 <- dat$p_het
  q <- ncol(B)
  nG <- n * G

  # If z is a vector, repeat it across components.
  # This is useful for the Gaussian case where z = y.
  if (!is.matrix(z)) {
    z <- matrix(as.numeric(z), nrow = n, ncol = G)
  } else {
    z <- as.matrix(z)
  }

  R1_list <- vector("list", G)
  E_list <- vector("list", G)
  b_list <- vector("list", G)

  Bhat <- matrix(0, nrow = nG, ncol = q)
  zhat <- numeric(nG)

  idx_g <- function(g) {
    ((g - 1) * n + 1):(g * n)
  }

  for (g in seq_len(G)) {
    ii <- idx_g(g)

    zg <- z[, g]
    wg <- pmax(w[, g], 0)
    sw <- sqrt(wg)

    Aw <- A * sw
    zw <- zg * sw

    qrA <- qr(Aw)

    R1g <- qr.R(qrA)[seq_len(p1), seq_len(p1), drop = FALSE]

    tmp_z <- qr.qty(qrA, zw)
    b_g <- tmp_z[seq_len(p1)]

    if (q > 0) {
      Bw <- B * sw

      tmp_B <- qr.qty(qrA, Bw)
      E_g <- tmp_B[seq_len(p1), , drop = FALSE]

      Bhat_g <- qr.resid(qrA, Bw)
      Bhat[ii, ] <- Bhat_g
    } else {
      E_g <- matrix(numeric(0), nrow = p1, ncol = 0)
    }

    zhat_g <- qr.resid(qrA, zw)
    zhat[ii] <- zhat_g

    R1_list[[g]] <- R1g
    E_list[[g]] <- E_g
    b_list[[g]] <- b_g
  }

  if (q > 0) {
    qrB <- qr(Bhat)
    R2 <- qr.R(qrB)[seq_len(q), seq_len(q), drop = FALSE]
    beta <- as.numeric(qr.coef(qrB, zhat)[seq_len(q)])
  } else {
    R2 <- matrix(numeric(0), nrow = 0, ncol = 0)
    beta <- numeric(0)
  }

  beta_g <- matrix(NA_real_, nrow = G, ncol = p1)

  for (g in seq_len(G)) {
    if (q > 0) {
      rhs <- b_list[[g]] - as.vector(E_list[[g]] %*% beta)
    } else {
      rhs <- b_list[[g]]
    }

    beta_g[g, ] <- backsolve(R1_list[[g]], rhs)
  }

  out <- list(
    beta_g = beta_g,
    beta = beta
  )


  if (isTRUE(control$return_qr_parts)) {
    out$R1_list <- R1_list
    out$E_list <- E_list
    out$R2 <- R2
  }

  out
}
