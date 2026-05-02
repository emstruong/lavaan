# lav_stan_data.R
#
# Convert a fitted-or-prepared lavaan model object into the data list
# consumed by inst/stan/lavaan_universal.stan.
#
# Strategy: walk lavmodel@GLIST and use lavmodel@m.free.idx /
# lavmodel@x.free.idx -- the existing sparse triplets that lavaan already
# maintains for its own x <-> matrix mapping. This guarantees we agree
# with lav_model_set_parameters() on parameter ordering.
#
# Scope (M1): single-group continuous data, no mean structure, LISREL.
# Returns NULL with a note if asked to convert a model outside scope.

# Find the matrix index in GLIST corresponding to a named LISREL matrix
# for a given group (1-based). Returns NA_integer_ if the matrix is not
# present in this group.
lav_stan_glist_idx <- function(lavmodel, group_idx, mat_name) {
  nmat       <- lavmodel@nmat
  start_idx  <- cumsum(c(0L, nmat))[group_idx] + 1L
  end_idx    <- start_idx + nmat[group_idx] - 1L
  full_names <- names(lavmodel@GLIST)[start_idx:end_idx]
  hit        <- which(full_names == mat_name)
  if (length(hit) == 0L) {
    return(NA_integer_)
  }
  start_idx + hit[1L] - 1L
}

# Extract triplets for a single LISREL matrix:
#   list(free_row, free_col, free_par,    # 1-based; par index into the
#                                         # full free-parameter vector x
#        fix_row,  fix_col,  fix_val)     # 1-based; fixed numeric values
lav_stan_matrix_triplets <- function(lavmodel, glist_index) {
  if (is.na(glist_index)) {
    return(list(
      free_row = integer(0), free_col = integer(0), free_par = integer(0),
      fix_row  = integer(0), fix_col  = integer(0), fix_val  = numeric(0)
    ))
  }
  mat        <- lavmodel@GLIST[[glist_index]]
  m_free_idx <- lavmodel@m.free.idx[[glist_index]]
  x_free_idx <- lavmodel@x.free.idx[[glist_index]]
  n_row      <- nrow(mat)

  # Free entries: (m_free_idx is the linear index into the matrix).
  if (length(m_free_idx) > 0L) {
    free_row <- ((m_free_idx - 1L) %% n_row) + 1L
    free_col <- ((m_free_idx - 1L) %/% n_row) + 1L
    free_par <- as.integer(x_free_idx)
  } else {
    free_row <- integer(0); free_col <- integer(0); free_par <- integer(0)
  }

  # Fixed entries: every nonzero, non-free element is "fixed". For a fully
  # zero matrix, this is empty.
  is_free <- logical(length(mat))
  if (length(m_free_idx) > 0L) is_free[m_free_idx] <- TRUE
  fix_lin <- which(!is_free & mat != 0)
  if (length(fix_lin) > 0L) {
    fix_row <- ((fix_lin - 1L) %% n_row) + 1L
    fix_col <- ((fix_lin - 1L) %/% n_row) + 1L
    fix_val <- as.numeric(mat[fix_lin])
  } else {
    fix_row <- integer(0); fix_col <- integer(0); fix_val <- numeric(0)
  }

  list(free_row = free_row, free_col = free_col, free_par = free_par,
       fix_row  = fix_row,  fix_col  = fix_col,  fix_val  = fix_val)
}

# Top-level conversion: lavaan objects -> Stan data list.
# Returns NULL if outside M1 scope (multigroup, mean structure, FIML, etc.).
lav_stan_data <- function(lavmodel = NULL,
                          lavsamplestats = NULL,
                          lavdata = NULL,
                          lavoptions = NULL) {
  # ---- M1 scope guards ----
  if (lavmodel@ngroups > 1L) {
    return(structure(NULL, reason =
      gettext("Multigroup models are not supported in this Stan backend yet.")))
  }
  if (lavmodel@meanstructure) {
    return(structure(NULL, reason =
      gettext("Mean structure is not supported in this Stan backend yet.")))
  }
  if (lavmodel@categorical) {
    return(structure(NULL, reason =
      gettext("Categorical data is not supported in this Stan backend yet.")))
  }
  if (!identical(lavmodel@representation, "LISREL")) {
    return(structure(NULL, reason =
      gettext("Only LISREL representation is supported in this Stan backend.")))
  }
  if (lavdata@nlevels > 1L) {
    return(structure(NULL, reason =
      gettext("Multilevel models are not supported in this Stan backend yet.")))
  }
  miss_arg <- lavoptions$missing
  if (!is.null(miss_arg) && !miss_arg %in% c("listwise", "default", "none")) {
    return(structure(NULL, reason =
      gettext("FIML is not supported in this Stan backend yet.")))
  }
  if (lavmodel@composites) {
    return(structure(NULL, reason =
      gettext("Composite constructs are not supported in this Stan backend yet.")))
  }
  if (length(lavmodel@ceq.nonlinear.idx) > 0L ||
      length(lavmodel@cin.nonlinear.idx) > 0L) {
    return(structure(NULL, reason =
      gettext("Nonlinear constraints are not supported in this Stan backend yet.")))
  }

  group_idx  <- 1L
  n_var      <- lavmodel@nvar[group_idx]
  lambda_idx <- lav_stan_glist_idx(lavmodel, group_idx, "lambda")
  theta_idx  <- lav_stan_glist_idx(lavmodel, group_idx, "theta")
  psi_idx    <- lav_stan_glist_idx(lavmodel, group_idx, "psi")
  beta_idx   <- lav_stan_glist_idx(lavmodel, group_idx, "beta")

  if (is.na(psi_idx)) {
    n_fac <- 0L
  } else {
    n_fac <- nrow(lavmodel@GLIST[[psi_idx]])
  }
  if (n_fac < 1L) {
    return(structure(NULL, reason =
      gettext("Pure path models (no latent factors) are not supported in this Stan backend yet.")))
  }

  lambda_t <- lav_stan_matrix_triplets(lavmodel, lambda_idx)
  theta_t  <- lav_stan_matrix_triplets(lavmodel, theta_idx)
  psi_t    <- lav_stan_matrix_triplets(lavmodel, psi_idx)
  beta_t   <- lav_stan_matrix_triplets(lavmodel, beta_idx)

  # Drop redundant upper-triangle entries from THETA / PSI (Stan code
  # symmetrises automatically). Keep diagonal.
  keep_lower <- function(t) {
    keep <- t$free_row >= t$free_col
    t$free_row <- t$free_row[keep]
    t$free_col <- t$free_col[keep]
    t$free_par <- t$free_par[keep]
    keep_f <- t$fix_row >= t$fix_col
    t$fix_row <- t$fix_row[keep_f]
    t$fix_col <- t$fix_col[keep_f]
    t$fix_val <- t$fix_val[keep_f]
    t
  }
  theta_t <- keep_lower(theta_t)
  psi_t   <- keep_lower(psi_t)

  use_beta <- if (is.na(beta_idx)) 0L else 1L

  list(
    n_var          = n_var,
    n_fac          = as.integer(n_fac),
    n_obs          = as.integer(lavsamplestats@nobs[[group_idx]]),
    n_par          = lavmodel@nx.free,
    sample_cov     = lavsamplestats@cov[[group_idx]],

    lambda_n_free  = length(lambda_t$free_par),
    lambda_free_row = as.array(as.integer(lambda_t$free_row)),
    lambda_free_col = as.array(as.integer(lambda_t$free_col)),
    lambda_free_par = as.array(as.integer(lambda_t$free_par)),
    lambda_n_fix   = length(lambda_t$fix_val),
    lambda_fix_row = as.array(as.integer(lambda_t$fix_row)),
    lambda_fix_col = as.array(as.integer(lambda_t$fix_col)),
    lambda_fix_val = as.array(as.numeric(lambda_t$fix_val)),

    theta_n_free   = length(theta_t$free_par),
    theta_free_row = as.array(as.integer(theta_t$free_row)),
    theta_free_col = as.array(as.integer(theta_t$free_col)),
    theta_free_par = as.array(as.integer(theta_t$free_par)),
    theta_n_fix    = length(theta_t$fix_val),
    theta_fix_row  = as.array(as.integer(theta_t$fix_row)),
    theta_fix_col  = as.array(as.integer(theta_t$fix_col)),
    theta_fix_val  = as.array(as.numeric(theta_t$fix_val)),

    psi_n_free     = length(psi_t$free_par),
    psi_free_row   = as.array(as.integer(psi_t$free_row)),
    psi_free_col   = as.array(as.integer(psi_t$free_col)),
    psi_free_par   = as.array(as.integer(psi_t$free_par)),
    psi_n_fix      = length(psi_t$fix_val),
    psi_fix_row    = as.array(as.integer(psi_t$fix_row)),
    psi_fix_col    = as.array(as.integer(psi_t$fix_col)),
    psi_fix_val    = as.array(as.numeric(psi_t$fix_val)),

    use_beta       = use_beta,
    beta_n_free    = length(beta_t$free_par),
    beta_free_row  = as.array(as.integer(beta_t$free_row)),
    beta_free_col  = as.array(as.integer(beta_t$free_col)),
    beta_free_par  = as.array(as.integer(beta_t$free_par)),
    beta_n_fix     = length(beta_t$fix_val),
    beta_fix_row   = as.array(as.integer(beta_t$fix_row)),
    beta_fix_col   = as.array(as.integer(beta_t$fix_col)),
    beta_fix_val   = as.array(as.numeric(beta_t$fix_val))
  )
}
