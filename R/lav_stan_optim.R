# lav_stan_optim.R
#
# Driver for optim.method = "stan". Returns the same attribute-rich numeric
# vector that lav_model_estimate() returns, so step12+ downstream code is
# unchanged.

# Backend-agnostic optimize call: returns list(par, fx, return_code, raw)
lav_stan_optimize <- function(compiled,
                              data_list,
                              init = NULL,
                              jacobian = FALSE,
                              algorithm = "lbfgs",
                              quiet = TRUE) {
  if (identical(compiled$backend, "cmdstanr")) {
    init_arg <- if (is.null(init)) NULL else list(theta_par = init)
    fit <- compiled$model$optimize(
      data      = data_list,
      init      = init_arg,
      jacobian  = jacobian,
      algorithm = algorithm,
      refresh   = if (quiet) 0 else 100
    )
    par_df <- fit$mle()
    par_vec <- as.numeric(par_df[grep("^theta_par", names(par_df))])
    fx_val  <- as.numeric(fit$lp())
    rc      <- as.integer(fit$return_codes())
    return(list(par = par_vec, fx = fx_val, return_code = rc, raw = fit))
  }
  if (identical(compiled$backend, "rstan")) {
    init_arg <- if (is.null(init)) "random" else list(theta_par = init)
    fit <- rstan::optimizing(
      object  = compiled$model,
      data    = data_list,
      init    = init_arg,
      hessian = FALSE,
      verbose = !quiet,
      algorithm = toupper(algorithm)
    )
    par_vec <- as.numeric(fit$par[grep("^theta_par", names(fit$par))])
    return(list(par = par_vec, fx = as.numeric(fit$value),
                return_code = as.integer(fit$return_code), raw = fit))
  }
  lav_msg_fixme("Unknown backend in lav_stan_optimize: ", compiled$backend)
}

# Top-level entry point invoked from step11 when optim.method = "stan".
# Mirrors the contract of lav_model_estimate(): returns a numeric vector x
# with attributes (iterations, converged, warn.txt, control, dx, fx).
lav_optim_stan <- function(lavmodel = NULL,
                           lavpartable = NULL,
                           lavsamplestats = NULL,
                           lavdata = NULL,
                           lavoptions = NULL,
                           lavcache = NULL) {
  backend <- lav_stan_check_install(
    backend = lavoptions$stan.backend %||% "auto", quiet = FALSE
  )
  if (is.null(backend)) {
    lav_msg_stop(gettext(
      "Stan backend selected via optim.method= \"stan\" but no Stan
      installation was found. Install cmdstanr and run
      cmdstanr::install_cmdstan(), or set optim.method= \"nlminb\"."
    ))
  }

  data_list <- lav_stan_data(lavmodel = lavmodel,
                             lavsamplestats = lavsamplestats,
                             lavdata = lavdata,
                             lavoptions = lavoptions)
  if (is.null(data_list)) {
    reason <- attr(data_list, "reason", exact = TRUE)
    lav_msg_stop(gettextf(
      "optim.method= \"stan\" cannot be used for this model: %1$s",
      reason
    ))
  }

  compiled <- lav_stan_compile(backend = backend, quiet = TRUE)

  # Starting values: lavaan already computed them and stored them in
  # lavmodel via @x.start; pull them out via partable$start.
  start_x <- lavpartable$start[lavpartable$free > 0L]
  start_x <- start_x[order(lavpartable$free[lavpartable$free > 0L])]

  res <- lav_stan_optimize(
    compiled  = compiled,
    data_list = data_list,
    init      = as.numeric(start_x),
    jacobian  = FALSE,
    algorithm = "lbfgs",
    quiet     = !lav_verbose()
  )

  x <- res$par

  # Re-evaluate the lavaan objective at the converged x so that step12+
  # see the same fx convention (e.g. -2*logL/N for ML) that the nlminb
  # path produces. This keeps fit measures and test stats correct
  # regardless of how Stan parameterises its target.
  glist_at_x <- lav_model_x2glist(lavmodel, x = x)
  fx <- lav_model_objective(
    lavmodel       = lavmodel,
    GLIST          = glist_at_x,
    lavsamplestats = lavsamplestats,
    lavdata        = lavdata,
    lavcache       = lavcache
  )

  attr(x, "iterations") <- NA_integer_
  attr(x, "converged")  <- isTRUE(res$return_code == 0L)
  attr(x, "warn.txt")   <- ""
  attr(x, "control")    <- list(stan.backend = backend)
  attr(x, "dx")         <- rep(NA_real_, length(x))
  attr(x, "fx")         <- fx
  attr(x, "stan.fit")   <- res$raw

  x
}

# Local null-coalesce; lavaan does not import rlang.
`%||%` <- function(a, b) if (is.null(a)) b else a
