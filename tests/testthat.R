# testthat entry point. Discovered automatically by R CMD check.
if (requireNamespace("testthat", quietly = TRUE)) {
  library(testthat)
  library(lavaan)
  test_check("lavaan")
}
