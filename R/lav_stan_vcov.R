# lav_stan_vcov.R
#
# Standard-error / vcov computation for fits produced by the Stan backend.
#
# M1: returns NULL so callers fall back to lavaan's existing handling
# (typically: vcov filled with NA, summary still prints estimates). The
# user must use se = "none" to get a clean summary in M1.
#
# Roadmap:
#   M4: $hessian()-based observed information; $laplace() draws for CIs
#   M5: $expose_model_methods() + $grad_log_prob() for per-case scores,
#       feeding lavaan's existing sandwich code in lav_model_vcov.R so
#       MLR robust SEs continue to work bit-for-bit.

lav_stan_vcov <- function(lavmodel = NULL,
                          lavoptions = NULL,
                          stan_fit = NULL) {
  if (is.null(stan_fit)) return(NULL)
  # Placeholder for M4-M5.
  NULL
}
