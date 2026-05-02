// lavaan universal Stan program
//
// Scope (M1): single-group, continuous indicators, complete data,
// no mean structure, LISREL representation, ML estimator, n_fac >= 1.
//
// Roadmap:
//   M1.5: pure path models (n_fac = 0)
//   M2:   mean structure (NU, ALPHA assembly + (ybar - mu) term)
//   M3:   multigroup (outer group loop) and FIML (per-pattern slicing)
//   M6:   composite-construct residual variance assembly
//
// Wire protocol:
//   The R side builds sparse "triplet" arrays mapping free parameter indices
//   (1..n_par) to (row, col) slots of LAMBDA / THETA / PSI / BETA, plus
//   arrays of fixed (row, col, value) entries. The Stan side scatter-
//   assembles the matrices each iteration. This avoids any per-model
//   code generation, so the Stan binary is compiled once and reused
//   across every lavaan model.

functions {
  // Sufficient-statistic multivariate normal log-likelihood:
  //   -0.5 * n * ( p*log(2pi) + log|Sigma| + tr(Sigma^{-1} * S) )
  // where S is the sample covariance computed with N (lavaan default).
  real lav_mvn_suffstat_lpdf(matrix sample_cov, int n_obs,
                             matrix sigma_hat) {
    int p = rows(sigma_hat);
    matrix[p, p] L = cholesky_decompose(sigma_hat);
    real log_det = 2 * sum(log(diagonal(L)));
    matrix[p, p] sigma_inv_S = mdivide_left_spd(sigma_hat, sample_cov);
    real tr_term = trace(sigma_inv_S);
    return -0.5 * n_obs * (p * log(2 * pi()) + log_det + tr_term);
  }
}

data {
  // ---- dimensions ----
  int<lower=1> n_var;          // number of observed variables
  int<lower=1> n_fac;          // number of latent factors (>= 1 in M1)
  int<lower=1> n_obs;          // sample size
  int<lower=1> n_par;          // number of free parameters

  // ---- sample sufficient statistics ----
  matrix[n_var, n_var] sample_cov;

  // ---- LAMBDA (n_var x n_fac) ----
  int<lower=0> lambda_n_free;
  array[lambda_n_free] int<lower=1> lambda_free_row;
  array[lambda_n_free] int<lower=1> lambda_free_col;
  array[lambda_n_free] int<lower=1> lambda_free_par;

  int<lower=0> lambda_n_fix;
  array[lambda_n_fix] int<lower=1> lambda_fix_row;
  array[lambda_n_fix] int<lower=1> lambda_fix_col;
  vector[lambda_n_fix] lambda_fix_val;

  // ---- THETA (n_var x n_var, symmetric, R sends lower triangle) ----
  int<lower=0> theta_n_free;
  array[theta_n_free] int<lower=1> theta_free_row;
  array[theta_n_free] int<lower=1> theta_free_col;
  array[theta_n_free] int<lower=1> theta_free_par;

  int<lower=0> theta_n_fix;
  array[theta_n_fix] int<lower=1> theta_fix_row;
  array[theta_n_fix] int<lower=1> theta_fix_col;
  vector[theta_n_fix] theta_fix_val;

  // ---- PSI (n_fac x n_fac, symmetric, R sends lower triangle) ----
  int<lower=0> psi_n_free;
  array[psi_n_free] int<lower=1> psi_free_row;
  array[psi_n_free] int<lower=1> psi_free_col;
  array[psi_n_free] int<lower=1> psi_free_par;

  int<lower=0> psi_n_fix;
  array[psi_n_fix] int<lower=1> psi_fix_row;
  array[psi_n_fix] int<lower=1> psi_fix_col;
  vector[psi_n_fix] psi_fix_val;

  // ---- BETA (n_fac x n_fac, structural) ----
  int<lower=0, upper=1> use_beta;
  int<lower=0> beta_n_free;
  array[beta_n_free] int<lower=1> beta_free_row;
  array[beta_n_free] int<lower=1> beta_free_col;
  array[beta_n_free] int<lower=1> beta_free_par;

  int<lower=0> beta_n_fix;
  array[beta_n_fix] int<lower=1> beta_fix_row;
  array[beta_n_fix] int<lower=1> beta_fix_col;
  vector[beta_n_fix] beta_fix_val;
}

parameters {
  vector[n_par] theta_par;
}

transformed parameters {
  matrix[n_var, n_fac] lambda_mat = rep_matrix(0, n_var, n_fac);
  matrix[n_var, n_var] theta_mat  = rep_matrix(0, n_var, n_var);
  matrix[n_fac, n_fac] psi_mat    = rep_matrix(0, n_fac, n_fac);
  matrix[n_var, n_var] sigma_hat;

  // Scatter-assemble LAMBDA
  for (i in 1:lambda_n_fix) {
    lambda_mat[lambda_fix_row[i], lambda_fix_col[i]] = lambda_fix_val[i];
  }
  for (i in 1:lambda_n_free) {
    lambda_mat[lambda_free_row[i], lambda_free_col[i]]
      = theta_par[lambda_free_par[i]];
  }

  // Scatter-assemble THETA (symmetric)
  for (i in 1:theta_n_fix) {
    theta_mat[theta_fix_row[i], theta_fix_col[i]] = theta_fix_val[i];
    if (theta_fix_row[i] != theta_fix_col[i]) {
      theta_mat[theta_fix_col[i], theta_fix_row[i]] = theta_fix_val[i];
    }
  }
  for (i in 1:theta_n_free) {
    theta_mat[theta_free_row[i], theta_free_col[i]]
      = theta_par[theta_free_par[i]];
    if (theta_free_row[i] != theta_free_col[i]) {
      theta_mat[theta_free_col[i], theta_free_row[i]]
        = theta_par[theta_free_par[i]];
    }
  }

  // Scatter-assemble PSI (symmetric)
  for (i in 1:psi_n_fix) {
    psi_mat[psi_fix_row[i], psi_fix_col[i]] = psi_fix_val[i];
    if (psi_fix_row[i] != psi_fix_col[i]) {
      psi_mat[psi_fix_col[i], psi_fix_row[i]] = psi_fix_val[i];
    }
  }
  for (i in 1:psi_n_free) {
    psi_mat[psi_free_row[i], psi_free_col[i]]
      = theta_par[psi_free_par[i]];
    if (psi_free_row[i] != psi_free_col[i]) {
      psi_mat[psi_free_col[i], psi_free_row[i]]
        = theta_par[psi_free_par[i]];
    }
  }

  // Sigma_hat = LAMBDA * (I - BETA)^{-1} * PSI * (I - BETA)^{-T} * LAMBDA' + THETA
  if (use_beta == 1) {
    matrix[n_fac, n_fac] beta_mat = rep_matrix(0, n_fac, n_fac);
    for (i in 1:beta_n_fix) {
      beta_mat[beta_fix_row[i], beta_fix_col[i]] = beta_fix_val[i];
    }
    for (i in 1:beta_n_free) {
      beta_mat[beta_free_row[i], beta_free_col[i]]
        = theta_par[beta_free_par[i]];
    }
    matrix[n_fac, n_fac] ib = diag_matrix(rep_vector(1.0, n_fac)) - beta_mat;
    matrix[n_fac, n_fac] ib_inv = inverse(ib);
    matrix[n_fac, n_fac] eta_cov = ib_inv * psi_mat * ib_inv';
    sigma_hat = lambda_mat * eta_cov * lambda_mat' + theta_mat;
  } else {
    sigma_hat = lambda_mat * psi_mat * lambda_mat' + theta_mat;
  }
}

model {
  target += lav_mvn_suffstat_lpdf(sample_cov | n_obs, sigma_hat);
}
