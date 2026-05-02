# Skip a test when no Stan backend is reachable. We treat "no backend" as
# a normal CRAN/CI condition rather than a test failure: lavaan must work
# without Stan installed.
skip_if_no_stan <- function(backend = "auto") {
  resolved <- try(
    lav_stan_backend_resolve(backend),
    silent = TRUE
  )
  if (inherits(resolved, "try-error")) {
    testthat::skip("No Stan backend (cmdstanr/rstan) available.")
  }
  invisible(resolved)
}

# Common CFA fixture used across regression tests.
hs_cfa_model <- "
  visual  =~ x1 + x2 + x3
  textual =~ x4 + x5 + x6
  speed   =~ x7 + x8 + x9
"
