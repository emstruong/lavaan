# Regression tests for the optim.method = "stan" backend.
#
# These tests run only when a Stan backend is reachable; otherwise they
# skip silently. The intent: when Stan IS available, we lock in that the
# Stan backend produces results indistinguishable (within numerical
# tolerance) from lavaan's default nlminb path on a small, fixed set of
# fixtures.
#
# When you add a milestone (mean structure, multigroup, FIML, composites),
# add a corresponding test_that() block here that exercises the new
# capability AGAINST the nlminb baseline.

test_that("backend resolver picks the right thing", {
  # No-Stan path should always succeed conceptually:
  # `lav_stan_backend_available("cmdstanr")` and `lav_stan_backend_available("rstan")`
  # both return logicals without erroring.
  expect_type(lav_stan_backend_available("cmdstanr"), "logical")
  expect_type(lav_stan_backend_available("rstan"),    "logical")

  # Bogus backend name -> error.
  expect_error(
    lav_stan_backend_resolve("not-a-real-backend"),
    regexp = "stan.backend"
  )
})

test_that("M1 scope guards refuse out-of-scope models gracefully", {
  # Multigroup must be refused (M3 territory).
  fit_baseline <- cfa(hs_cfa_model, data = HolzingerSwineford1939,
                     group = "school")
  data_list <- lav_stan_data(
    lavmodel       = fit_baseline@Model,
    lavsamplestats = fit_baseline@SampleStats,
    lavdata        = fit_baseline@Data,
    lavoptions     = fit_baseline@Options
  )
  expect_null(data_list)
  expect_match(attr(data_list, "reason"), regexp = "[Mm]ultigroup")
})

test_that("Stan ML matches nlminb ML on HS single-group CFA", {
  skip_if_no_stan()

  fit_nlminb <- cfa(hs_cfa_model, data = HolzingerSwineford1939,
                    se = "none")
  fit_stan   <- cfa(hs_cfa_model, data = HolzingerSwineford1939,
                    se = "none", optim.method = "stan")

  # Coefficients should agree to 4 decimals (Stan L-BFGS terminates on
  # `tol_rel_obj`; lavaan's nlminb on parameter and objective tolerances).
  expect_equal(coef(fit_stan), coef(fit_nlminb), tolerance = 1e-4)

  # Log-likelihood should agree to 6 decimals (the objective is the same
  # function in both cases; only the optimizer differs).
  expect_equal(as.numeric(logLik(fit_stan)),
               as.numeric(logLik(fit_nlminb)),
               tolerance = 1e-6)
})

test_that("Stan ML matches nlminb ML on a structural CFA-with-paths", {
  skip_if_no_stan()
  # Two-factor model with a latent regression -- exercises the BETA branch
  # of the Stan program.
  struct_model <- "
    visual  =~ x1 + x2 + x3
    textual =~ x4 + x5 + x6
    textual ~ visual
  "
  fit_nlminb <- sem(struct_model, data = HolzingerSwineford1939, se = "none")
  fit_stan   <- sem(struct_model, data = HolzingerSwineford1939, se = "none",
                    optim.method = "stan")
  expect_equal(coef(fit_stan), coef(fit_nlminb), tolerance = 1e-4)
})
