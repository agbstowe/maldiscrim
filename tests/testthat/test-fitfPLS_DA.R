# Load package built-in data required for testing
data("spectra100", package = "maldiscrim")

test_that("fit.fPLS_DA works correctly using B-splines", {
  # 1. Run the model (lowering n_folds to speed up the test execution)
  model <- fit.fPLS_DA(data = spectra100, method = "bsplines", nbasis = 50, n_folds = 2)

  # 2. Check expected structural outputs
  expect_s3_class(model, "fPLS_DA")
  expect_type(model, "list")
  expect_named(model, c("pls_model", "lda_model", "ncompOpt", "const_idx", "method", "levels"))

  # 3. Check baseline mathematical/logical values
  expect_true(model$ncompOpt > 0)
  expect_equal(model$method, "bsplines")
})

test_that("fit.fPLS_DA throws appropriate errors for invalid parameters", {
  # Verify that the function gracefully halts when encountering critical input issues
  expect_error(fit.fPLS_DA(data = NULL))
  expect_error(fit.fPLS_DA(data = spectra100, method = "non_existent_method"))
  expect_error(fit.fPLS_DA(data = spectra100, n_folds = -5))
})
