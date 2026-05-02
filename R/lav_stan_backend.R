# lav_stan_backend.R
#
# Backend abstraction for the Stan-based optimizer.
#
# Two backends are anticipated:
#   * "cmdstanr" -- shipped first; CmdStan installed by the user via
#                   cmdstanr::install_cmdstan(). No install-time C++.
#   * "rstan"    -- planned (M7); runtime stan_model() compile, no
#                   rstantools precompile so we never add LinkingTo.
#
# Selection precedence:
#   1. lavoptions$stan.backend (per-call argument, "auto"/"cmdstanr"/"rstan")
#   2. getOption("lavaan.stan.backend", default = "auto")
#   3. auto -> first available among {cmdstanr, rstan}
#
# All Stan-touching code in lavaan calls through this layer; no other file
# should reference cmdstanr::/rstan:: directly.

# Returns one of "cmdstanr" or "rstan", or stops with an actionable error.
lav_stan_backend_resolve <- function(requested = "auto") {
  if (is.null(requested) || identical(requested, "auto")) {
    requested <- getOption("lavaan.stan.backend", default = "auto")
  }
  if (identical(requested, "auto")) {
    if (lav_stan_backend_available("cmdstanr")) return("cmdstanr")
    if (lav_stan_backend_available("rstan"))    return("rstan")
    lav_msg_stop(gettext(
      "Stan backend requested but neither cmdstanr nor rstan is installed."
    ))
  }
  if (!requested %in% c("cmdstanr", "rstan")) {
    lav_msg_stop(gettextf(
      "stan.backend= argument must be one of %1$s; got %2$s.",
      "\"auto\", \"cmdstanr\", \"rstan\"", dQuote(requested)
    ))
  }
  if (!lav_stan_backend_available(requested)) {
    lav_msg_stop(gettextf(
      "stan.backend= %1$s requested but the package is not installed.",
      dQuote(requested)
    ))
  }
  requested
}

# Returns TRUE if the named backend's R package is loadable AND its underlying
# Stan toolchain is reachable. For cmdstanr this means a CmdStan install is
# present; for rstan it means rstan can be attached.
lav_stan_backend_available <- function(backend) {
  if (identical(backend, "cmdstanr")) {
    if (!requireNamespace("cmdstanr", quietly = TRUE)) return(FALSE)
    cv <- try(cmdstanr::cmdstan_version(error_on_NA = FALSE), silent = TRUE)
    return(!inherits(cv, "try-error") && !is.na(cv))
  }
  if (identical(backend, "rstan")) {
    return(requireNamespace("rstan", quietly = TRUE))
  }
  FALSE
}

# Single point of contact for "is the Stan path usable at all?". Used by
# step11 to decide whether to fall back to nlminb with a note instead of
# erroring out hard.
lav_stan_check_install <- function(backend = "auto", quiet = FALSE) {
  resolved <- try(lav_stan_backend_resolve(backend), silent = TRUE)
  if (inherits(resolved, "try-error")) {
    if (!quiet) {
      lav_msg_warn(gettextf(
        "Stan backend not available: %1$s",
        conditionMessage(attr(resolved, "condition"))
      ))
    }
    return(NULL)
  }
  resolved
}
