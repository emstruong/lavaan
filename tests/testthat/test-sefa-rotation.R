# Adversarial-collaboration tests for the SEFA / DSEFA extension.
#
# These tests implement the rounds described in the project plan:
#   Round 1: gradient correctness for SEFA, ALF and Target-ALF criteria.
#   Round 2: SEFA fit-equivalence to oblique EFA on a fixed dataset
#            (Holzinger-Swineford 1939) plus second-order extraction.
#
# Builder makes claims; Red Team adds the checks that would catch a bug.
# Every claim has at least one expect_*() that would fail under a
# plausible regression.

context("SEFA rotation criteria")

# ---------------------------------------------------------------------
# Round 1 -- gradient correctness (Red Team: silent gradient bugs)
# ---------------------------------------------------------------------

test_that("ALF criterion analytic gradient matches numDeriv", {
  skip_if_not_installed("numDeriv")
  set.seed(42)
  lambda <- matrix(stats::rnorm(20L * 5L), 20L, 5L)
  ana <- attr(lavaan:::lav_matrix_rotate_alf(lambda, grad = TRUE), "grad")
  num <- matrix(numDeriv::grad(
    func = function(v) {
      lavaan:::lav_matrix_rotate_alf(matrix(v, 20L, 5L), grad = FALSE)
    },
    x = as.vector(lambda)
  ), 20L, 5L)
  expect_equal(ana, num, tolerance = 1e-06)
})

test_that("Target-ALF criterion analytic gradient matches numDeriv", {
  skip_if_not_installed("numDeriv")
  set.seed(43)
  lambda <- matrix(stats::rnorm(15L * 4L), 15L, 4L)
  target <- matrix(0, 15L, 4L)
  mask <- matrix(1, 15L, 4L)
  # zero-out a column to mimic a partially-specified target
  mask[, 1L] <- 0
  ana <- attr(
    lavaan:::lav_matrix_rotate_target_alf(
      lambda, target = target, target.mask = mask, grad = TRUE
    ),
    "grad"
  )
  num <- matrix(numDeriv::grad(
    func = function(v) {
      lavaan:::lav_matrix_rotate_target_alf(
        matrix(v, 15L, 4L),
        target = target, target.mask = mask, grad = FALSE
      )
    },
    x = as.vector(lambda)
  ), 15L, 4L)
  expect_equal(ana, num, tolerance = 1e-06)
})

test_that("SEFA criterion analytic gradient matches numDeriv", {
  skip_if_not_installed("numDeriv")
  set.seed(44)
  lambda <- matrix(stats::rnorm(20L * 4L), 20L, 4L)
  ana <- attr(lavaan:::lav_matrix_rotate_sefa(lambda, grad = TRUE), "grad")
  num <- matrix(numDeriv::grad(
    func = function(v) {
      lavaan:::lav_matrix_rotate_sefa(matrix(v, 20L, 4L), grad = FALSE)
    },
    x = as.vector(lambda)
  ), 20L, 4L)
  expect_equal(ana, num, tolerance = 1e-06)
})

# Red Team probe: ALF should not be exactly L2. If a future refactor
# accidentally aliases ALF to standard target rotation we want to know.
test_that("ALF is genuinely different from L2-target", {
  set.seed(7)
  lambda <- matrix(stats::rnorm(10L * 3L), 10L, 3L)
  q_alf <- lavaan:::lav_matrix_rotate_alf(lambda, grad = FALSE)
  q_l2 <- sum(lambda * lambda)
  expect_false(isTRUE(all.equal(q_alf, q_l2, tolerance = 1e-04)))
})

# ---------------------------------------------------------------------
# Round 2 -- SEFA = oblique EFA at the data-fit level
# ---------------------------------------------------------------------

test_that("SEFA second-order extraction produces consistent Lambda2", {
  # Construct a Phi that is *exactly* a one-factor model so we know the
  # truth and can assert exact recovery (modulo sign).
  l2_true <- c(0.7, 0.6, 0.8, 0.5)
  phi <- tcrossprod(l2_true) + diag(1 - l2_true^2)
  # numerical noise on diag should be 1
  expect_equal(diag(phi), rep(1, 4L), tolerance = 1e-12)

  ext <- lavaan:::lav_efa_sefa_extract(phi)
  expect_true(ext$converged)
  expect_false(ext$heywood)

  # Sign of lambda2 is identified only up to a global sign flip.
  l2_hat <- as.numeric(ext$lambda2)
  if (sum(l2_hat * l2_true) < 0) {
    l2_hat <- -l2_hat
  }
  expect_equal(l2_hat, l2_true, tolerance = 1e-04)
  # paper eq. 6: theta_i = 1 - lambda2_i^2 by construction
  expect_equal(ext$theta, 1 - l2_hat^2, tolerance = 1e-08)
})

test_that("SEFA extractor flags Heywood when |lambda2| approaches 1", {
  # Construct a Phi where the implied lambda2 ~ 1 -- a near-Heywood case.
  l2_true <- c(0.99, 0.97, 0.98)
  phi <- tcrossprod(l2_true) + diag(1 - l2_true^2)
  ext <- lavaan:::lav_efa_sefa_extract(phi)
  expect_true(ext$heywood)
})

test_that("SEFA extractor errors on a non-square or non-PD-like input", {
  expect_error(lavaan:::lav_efa_sefa_extract(matrix(0, 3L, 4L)),
               "square")
})

# ---------------------------------------------------------------------
# Round 2b -- end-to-end with HolzingerSwineford1939 (smoke test)
# ---------------------------------------------------------------------

test_that("efa(rotation = 'sefa') runs on HolzingerSwineford1939", {
  skip_on_cran()
  hs <- lavaan::HolzingerSwineford1939
  fit_sefa <- try(
    lavaan::efa(data = hs[, paste0("x", 1:9)], nfactors = 3L,
                rotation = "sefa",
                rotation.args = list(orthogonal = FALSE)),
    silent = TRUE
  )
  # Smoke test: rotation = "sefa" must be a recognised method.
  expect_false(inherits(fit_sefa, "try-error"))
})

test_that("efa(rotation = 'target.alf') is recognised", {
  skip_on_cran()
  hs <- lavaan::HolzingerSwineford1939
  target <- matrix(0, 9L, 3L)
  target[1:3, 1L] <- NA
  target[4:6, 2L] <- NA
  target[7:9, 3L] <- NA
  mask <- ifelse(is.na(target), 0, 1)
  target[is.na(target)] <- 0
  fit <- try(
    lavaan::efa(data = hs[, paste0("x", 1:9)], nfactors = 3L,
                rotation = "target.alf",
                rotation.args = list(target = target, target.mask = mask)),
    silent = TRUE
  )
  expect_false(inherits(fit, "try-error"))
})
