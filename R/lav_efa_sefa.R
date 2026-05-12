# SEFA second-order extraction utilities
#
# Asparouhov, T. & Muthen, B. (2026). A unification of second-order
# and bi-factor exploratory factor analysis (SEFA-DSEFA paper).
#
# Section 3 of the paper establishes that SEFA with m1 first-order
# factors has the same data-fit as oblique EFA with m1 factors. The
# rotation that drives the first-order solution is governed by the
# Geomin penalty on Lambda1 (paper eq. 4 / eq. 11). The second-order
# structure is then read off the factor correlation matrix Phi.
#
# Given an oblique-EFA's Phi, this file fits a single-factor model
#
#     Phi = Lambda2 %*% t(Lambda2) + diag(theta)
#
# subject to the SEFA scaling constraint Var(F) = 1, equivalently
#
#     theta_i = 1 - lambda2_i^2     (paper eq. 6)
#
# which guarantees diag(Phi) = 1 and avoids the second-order Heywood
# cases noted on paper p. 4 / Table on p. 6 (the original Var(xi)=1
# scaling produced 54% non-convergence at N=300 vs. 0% under Var(F)=1).
#
# YR/EM 2026 -- initial version (Phase 1 of the SEFA extension)

# Internal: solve for Lambda2 (a vector of length m1) that minimises
# sum of squared off-diagonal residuals between Phi and
# (Lambda2 Lambda2' + diag(1 - Lambda2^2)) using a constrained-Newton
# step with bounds |lambda2_i| < 1 - heywood_tol.
#
# This is an iterative least-squares fit; for m1 == 1 there is a
# trivial closed form (lambda2 = 1) and for m1 == 2 we fall back to a
# direct solve. For m1 >= 3 we use a coordinate-descent inner loop on
# lambda2_i with the others held fixed.
lav_efa_sefa_fit_phi <- function(phi,                        # nolint
                                 max_iter = 200L,
                                 tol = 1e-08,
                                 heywood_tol = 1e-03) {
  if (!is.matrix(phi)) {
    lav_msg_stop(gettext("phi= argument must be a matrix."))
  }
  if (nrow(phi) != ncol(phi)) {
    lav_msg_stop(gettext("phi= argument must be a square matrix."))
  }

  m1 <- nrow(phi)

  if (m1 < 2L) {
    return(list(
      lambda2 = matrix(1, 1L, 1L),
      theta = 0,
      converged = TRUE,
      iter = 0L,
      heywood = FALSE,
      residual = matrix(0, 1L, 1L)
    ))
  }

  # numerical floor on |lambda2_i| to keep theta_i = 1 - lambda2_i^2 from
  # going to zero / negative during the inner solve. Distinct from the
  # heywood_tol used for *reporting* Heywood: the floor is just to keep
  # the numerics well-conditioned, while heywood_tol governs the user-
  # facing diagnostic.
  numeric_floor <- 1e-08

  # initial values: signed sqrt of average off-diagonal column sum,
  # bounded away from 1 to avoid Heywood at the start.
  off <- phi
  diag(off) <- 0
  init <- sign(rowSums(off)) * sqrt(pmin(0.9, abs(rowMeans(off))))
  init[init == 0] <- 0.3
  lambda2 <- as.numeric(init)

  # off-diagonal coordinate-descent loop
  converged <- FALSE
  iter <- 0L
  for (iter in seq_len(max_iter)) {
    lambda2_old <- lambda2
    for (i in seq_len(m1)) {
      # for each i, solve a 1D least-squares problem on lambda2_i
      # using off-diagonal entries Phi_{ij} ~= lambda2_i * lambda2_j.
      others <- lambda2[-i]
      rhs <- phi[i, -i]
      denom <- sum(others * others)
      if (denom < .Machine$double.eps) {
        next
      }
      new_val <- sum(others * rhs) / denom
      # tight numerical floor so theta stays positive
      max_abs <- 1 - numeric_floor
      if (abs(new_val) > max_abs) {
        new_val <- sign(new_val) * max_abs
      }
      lambda2[i] <- new_val
    }
    if (max(abs(lambda2 - lambda2_old)) < tol) {
      converged <- TRUE
      break
    }
  }

  theta <- 1 - lambda2 * lambda2
  # Report Heywood when any second-order communality is at or above
  # 1 - heywood_tol (equivalently, theta_i <= heywood_tol). See paper
  # p. 6 for the second-order Heywood phenomenon.
  heywood <- any(theta <= heywood_tol)

  implied <- tcrossprod(lambda2) + diag(theta)
  residual <- phi - implied

  list(
    lambda2 = matrix(lambda2, m1, 1L),
    theta = theta,
    converged = converged,
    iter = iter,
    heywood = heywood,
    residual = residual
  )
}

# Public extractor: take a fitted lavaan EFA object and produce the
# SEFA second-order parameters. Works on the *rotated* factor
# correlation matrix Phi.
#
# x         : a fitted lavaan EFA object (oblique rotation), or a
#             plain p x p factor correlation matrix.
# block     : EFA block index when x is a fitted lavaan object with
#             multiple EFA blocks (default 1L).
# max_iter, tol, heywood_tol : passed through to lav_efa_sefa_fit_phi.
#
# Returns a list with class "lav_sefa":
#   lambda1   : first-order rotated loadings (Lambda1)
#   phi       : factor correlation matrix used for extraction
#   lambda2   : second-order loadings (m1 x 1)
#   theta     : residual variances of first-order factors
#               (= 1 - lambda2^2 by construction; paper eq. 6)
#   residual  : phi - (lambda2 lambda2' + diag(theta))
#   converged : logical, did the inner solver converge?
#   heywood   : logical, did any 1 - lambda2_i^2 fall below
#               heywood_tol? (See paper p. 6 for context.)
lav_efa_sefa_extract <- function(x, block = 1L,                  # nolint
                                 max_iter = 200L,
                                 tol = 1e-08,
                                 heywood_tol = 1e-03) {
  if (is.matrix(x)) {
    phi <- x
    lambda1 <- NULL
  } else if (inherits(x, "efaList")) {
    # efaList is a list of fitted lavaan objects, one per nfactors. The
    # 'block' argument here selects which entry; if names like "nf3" are
    # present we accept either an integer index or that name string.
    target_obj <- if (is.character(block)) x[[block]] else x[[block]]
    if (is.null(target_obj) || !inherits(target_obj, "lavaan")) {
      lav_msg_stop(gettextf(
        "block= %s does not select a fitted lavaan object inside efaList.",
        dQuote(as.character(block))))
    }
    return(lav_efa_sefa_extract(
      x = target_obj, block = 1L,
      max_iter = max_iter, tol = tol, heywood_tol = heywood_tol
    ))
  } else if (inherits(x, "lavaan")) {
    est <- lavInspect(x, "est")
    # est can be either a single list (one block) or a list of lists
    # (multiple groups / blocks). Normalize:
    if (!is.null(est$psi) && !is.null(est$lambda)) {
      psi <- est$psi
      lambda1 <- est$lambda
    } else {
      blk <- est[[block]]
      psi <- blk$psi
      lambda1 <- blk$lambda
    }
    # For an oblique EFA, the factor correlation matrix is the
    # standardized psi: phi_ij = psi_ij / sqrt(psi_ii psi_jj).
    sds <- sqrt(diag(psi))
    if (any(sds <= 0)) {
      lav_msg_stop(gettext(
        "Non-positive factor variance detected; cannot extract SEFA."))
    }
    phi <- psi / outer(sds, sds)
  } else {
    lav_msg_stop(gettextf(
      "x= argument must be a fitted lavaan object or a matrix; got %s.",
      dQuote(class(x)[1L])))
  }

  fit <- lav_efa_sefa_fit_phi(phi, max_iter = max_iter, tol = tol,
                              heywood_tol = heywood_tol)

  out <- list(
    lambda1 = lambda1,
    phi = phi,
    lambda2 = fit$lambda2,
    theta = fit$theta,
    residual = fit$residual,
    converged = fit$converged,
    iter = fit$iter,
    heywood = fit$heywood
  )
  class(out) <- c("lav_sefa", "list")
  out
}

# Print method: show a compact summary of the SEFA extraction.
print.lav_sefa <- function(x, digits = 3L, ...) {                # nolint
  cat("lavaan SEFA second-order extraction\n")
  cat("------------------------------------\n")
  cat(sprintf("Number of first-order factors: %d\n",
              length(x$theta)))
  cat(sprintf("Inner-solver converged:        %s (iter = %d)\n",
              x$converged, x$iter))
  if (x$heywood) {
    cat("Heywood case detected: at least one factor has |lambda2| ~ 1.\n")
  }
  cat("\nLambda2 (second-order loadings):\n")
  print(round(as.numeric(x$lambda2), digits))
  cat("\ntheta (residual variances of first-order factors):\n")
  print(round(x$theta, digits))
  cat("\nMax abs off-diagonal residual:",
      formatC(max(abs(x$residual - diag(diag(x$residual)))),
              format = "g", digits = digits), "\n")
  invisible(x)
}

# User-facing exported wrapper. Follows the lavaan style guide: exported
# functions other than the core fitters use 'lav' prefix + CamelCase.
lavSefaExtract <- function(object, block = 1L,                   # nolint
                           max.iter = 200L,                      # nolint
                           tol = 1e-08,
                           heywood.tol = 1e-03) {                # nolint
  lav_efa_sefa_extract(
    x = object,
    block = block,
    max_iter = max.iter,
    tol = tol,
    heywood_tol = heywood.tol
  )
}

