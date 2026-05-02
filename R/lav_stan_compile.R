# lav_stan_compile.R
#
# Compile / cache helpers for the universal Stan program.
#
# The .stan source ships in inst/stan/lavaan_universal.stan. We compile it
# once per (backend, Stan version, source SHA) and cache the resulting
# model object on disk under tools::R_user_dir("lavaan", "cache").
#
# Compilation is deferred until first use of the Stan backend so that
# users who never opt into optim.method = "stan" never pay any cost.

# Path to the bundled Stan source.
lav_stan_source_path <- function() {
  system.file("stan", "lavaan_universal.stan", package = "lavaan",
              mustWork = TRUE)
}

# A compiled-model handle. For cmdstanr this is a CmdStanModel; for rstan
# (future) it will be a stanmodel S4. Wrapped in a small list so the rest
# of the code is backend-agnostic.
lav_stan_compile <- function(backend = "auto",
                             stan_file = lav_stan_source_path(),
                             cache_dir = NULL,
                             quiet = TRUE) {
  backend <- lav_stan_backend_resolve(backend)
  if (is.null(cache_dir)) {
    cache_dir <- tools::R_user_dir("lavaan", which = "cache")
  }
  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  }
  if (identical(backend, "cmdstanr")) {
    return(lav_stan_compile_cmdstanr(stan_file, cache_dir, quiet))
  }
  if (identical(backend, "rstan")) {
    return(lav_stan_compile_rstan(stan_file, cache_dir, quiet))
  }
  lav_msg_fixme("Unreachable backend in lav_stan_compile: ", backend)
}

lav_stan_compile_cmdstanr <- function(stan_file, cache_dir, quiet) {
  # cmdstanr already caches by stanc hash in dir=; we just give it our cache.
  mod <- cmdstanr::cmdstan_model(
    stan_file = stan_file,
    dir       = cache_dir,
    quiet     = quiet,
    compile   = TRUE
  )
  list(backend = "cmdstanr", model = mod, source = stan_file)
}

lav_stan_compile_rstan <- function(stan_file, cache_dir, quiet) {
  # M7: runtime stan_model() compile, no rstantools precompile.
  # Cache the stanmodel object as an .rds keyed by source SHA + rstan version.
  src   <- readLines(stan_file)
  sha   <- substr(rlang_sha1_or_digest(paste(src, collapse = "\n")), 1, 12)
  rver  <- as.character(utils::packageVersion("rstan"))
  cfile <- file.path(cache_dir, paste0("rstan_", sha, "_", rver, ".rds"))
  if (file.exists(cfile)) {
    mod <- readRDS(cfile)
  } else {
    mod <- rstan::stan_model(file = stan_file, verbose = !quiet)
    saveRDS(mod, cfile)
  }
  list(backend = "rstan", model = mod, source = stan_file)
}

# Tiny shim: prefer digest::digest if available, else a poor-man's hash.
# (lavaan does not currently depend on digest; we don't want to add it just
# for cache keys.)
rlang_sha1_or_digest <- function(x) {
  if (requireNamespace("digest", quietly = TRUE)) {
    return(digest::digest(x, algo = "sha1"))
  }
  # Fallback: collapse to length-bounded string. Not cryptographic; only
  # used as a cache key against the bundled .stan source bytes.
  paste0("len", nchar(x), "_chk", sum(utf8ToInt(substr(x, 1, 4096))))
}
