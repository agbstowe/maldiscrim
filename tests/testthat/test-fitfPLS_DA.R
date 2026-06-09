# Load package built-in data required for testing
data("spectra100", package = "maldiscrim")

# Fixed to match EXACTLY the names defined in your fitFPLS_DA function output list
expected_names <- c("pls_model", "lda_model", "method", "basis", "wavelets","ncompOpt", "const_idx",
  "levels", "labels", "argvals", "filter", "nlevels", "level", "nbasis", "boundary")

# B-splines

test_that("fitFPLS_DA returns a valid FPLS_DA object with B-splines", {
  model <- fitFPLS_DA(data = spectra100, method = "bsplines", nbasis = 50, nfolds = 2)

  # Class and structure (Matches class(model_obj) <- "FPLS_DA")
  expect_s3_class(model, "FPLS_DA")
  expect_type(model, "list")
  expect_named(model, expected_names)

  expect_equal(model$method, "bsplines")

  expect_equal(model$nbasis, 50)

  expect_true(model$ncompOpt > 0)
  expect_true(length(model$levels) > 1)

  expect_true(is.factor(model$labels))
  expect_true(is.logical(model$const_idx))

  # B-splines fields present, wavelets fields absent
  expect_false(is.null(model$basis))
  expect_false(is.null(model$argvals))
  expect_null(model$wavelets)
  expect_null(model$filter)
  expect_null(model$nlevels)
  expect_null(model$level)
  expect_null(model$boundary)
})


# Precomputed shortcut

test_that("fitFPLS_DA returns a valid object with precomputed = TRUE", {
  model <- fitFPLS_DA(precomputed = TRUE)
  expect_s3_class(model, "FPLS_DA")
  expect_type(model, "list")
  expect_named(model, expected_names)
})


# Input validation (guards)

test_that("fitFPLS_DA throws appropriate errors for invalid inputs", {

  # data
  expect_error(fitFPLS_DA(data = NULL), "'data' must be a matrix or a data.frame.")
  expect_error(fitFPLS_DA(data = "not_a_matrix"),  "'data' must be a matrix or a data.frame.")

  # method
  expect_error(fitFPLS_DA(data = spectra100, method = "non_existent_method"))

  # nbasis
  expect_error(fitFPLS_DA(data = spectra100, method = "bsplines", nbasis = -1), "'nbasis' must be a single positive integer.")
  expect_error(fitFPLS_DA(data = spectra100, method = "bsplines", nbasis = 1.5), "'nbasis' must be a single positive integer.")
  expect_error(fitFPLS_DA(data = spectra100, method = "bsplines", nbasis = 0),  "'nbasis' must be a single positive integer.")

  # ncomp
  expect_error(fitFPLS_DA(data = spectra100, method = "bsplines", nbasis = 50, ncomp = -1),
               "'ncomp' must be a single positive integer or NULL.")
  expect_error(fitFPLS_DA(data = spectra100, method = "bsplines",  nbasis = 50, ncomp = 2.5),
               "'ncomp' must be a single positive integer or NULL.")

  # nfolds
  expect_error(fitFPLS_DA(data = spectra100, method = "bsplines", nbasis = 50, nfolds = -5),
               "'nfolds' must be a single integer >= 2.")
  expect_error(fitFPLS_DA(data = spectra100, method = "bsplines",  nbasis = 50, nfolds = 1),
               "'nfolds' must be a single integer >= 2.")
  expect_error(fitFPLS_DA(data = spectra100, method = "bsplines",  nbasis = 50, nfolds = 2.5),
               "'nfolds' must be a single integer >= 2.")
})
