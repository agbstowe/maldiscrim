# Load package built-in data required for testing
data("spectra100", package = "maldiscrim")

# Expected names of the FPLS_DA model object
expected_names <- c("pls_model", "lda_model", "method", "basis", "wavelets",
  "ncompOpt", "const_idx", "levels", "labels", "argvals","filter", "nlevels", "level", "nbasis", "boundary")


# B-splines

test_that("fitFPLS_DA returns a valid FPLS_DA object with B-splines", {
  model <- fitFPLS_DA(data = spectra100, method = "bsplines", nbasis = 50, nfolds = 2)

  # Class and structure
  expect_s3_class(model, "FPLS_DA")
  expect_type(model, "list")
  expect_named(model, expected_names)

  # Method coherence
  expect_equal(model$method, "bsplines")

  # Parameter coherence
  expect_equal(model$nbasis, 50)

  # Mathematical validity
  expect_true(model$ncompOpt > 0)
  expect_true(length(model$levels) > 1)

  # Type checks
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


# -- Wavelets -----------------------------------------------------------------

test_that("fitFPLS_DA returns a valid FPLS_DA object with wavelets", {
  model <- fitFPLS_DA(data = spectra100, method = "wavelets",
                      level = 3, nfolds = 2)

  # Class and structure
  expect_s3_class(model, "FPLS_DA")
  expect_type(model, "list")
  expect_named(model, expected_names)

  # Method coherence
  expect_equal(model$method, "wavelets")

  # Mathematical validity
  expect_true(model$ncompOpt > 0)
  expect_true(length(model$levels) > 1)

  # Type checks
  expect_true(is.factor(model$labels))
  expect_true(is.logical(model$const_idx))

  # Wavelets fields present, B-splines fields absent
  expect_false(is.null(model$wavelets))
  expect_false(is.null(model$filter))
  expect_false(is.null(model$nlevels))
  expect_false(is.null(model$level))
  expect_false(is.null(model$boundary))
  expect_null(model$basis)
  expect_null(model$argvals)
  expect_null(model$nbasis)
})


# -- Precomputed shortcut -----------------------------------------------------

test_that("fitFPLS_DA returns a valid object with precomputed = TRUE", {
  model <- fitFPLS_DA(precomputed = TRUE)
  expect_s3_class(model, "FPLS_DA")
  expect_type(model, "list")
  expect_named(model, expected_names)
})


# -- Input validation (guards) ------------------------------------------------

test_that("fitFPLS_DA throws appropriate errors for invalid inputs", {

  # data
  expect_error(fitFPLS_DA(data = NULL),
               "'data' must be a matrix or a data.frame.")
  expect_error(fitFPLS_DA(data = "not_a_matrix"),
               "'data' must be a matrix or a data.frame.")

  # method
  expect_error(fitFPLS_DA(data = spectra100, method = "non_existent_method"))

  # nbasis
  expect_error(fitFPLS_DA(data = spectra100, method = "bsplines", nbasis = -1),
               "'nbasis' must be a single positive integer.")
  expect_error(fitFPLS_DA(data = spectra100, method = "bsplines", nbasis = 1.5),
               "'nbasis' must be a single positive integer.")
  expect_error(fitFPLS_DA(data = spectra100, method = "bsplines", nbasis = 0),
               "'nbasis' must be a single positive integer.")

  # ncomp
  expect_error(fitFPLS_DA(data = spectra100, method = "bsplines",
                          nbasis = 50, ncomp = -1),
               "'ncomp' must be a single positive integer or NULL.")
  expect_error(fitFPLS_DA(data = spectra100, method = "bsplines",
                          nbasis = 50, ncomp = 2.5),
               "'ncomp' must be a single positive integer or NULL.")

  # nfolds
  expect_error(fitFPLS_DA(data = spectra100, method = "bsplines",
                          nbasis = 50, nfolds = -5),
               "'nfolds' must be a single integer >= 2.")
  expect_error(fitFPLS_DA(data = spectra100, method = "bsplines",
                          nbasis = 50, nfolds = 1),
               "'nfolds' must be a single integer >= 2.")
  expect_error(fitFPLS_DA(data = spectra100, method = "bsplines",
                          nbasis = 50, nfolds = 2.5),
               "'nfolds' must be a single integer >= 2.")

  # nlevels
  expect_error(fitFPLS_DA(data = spectra100, method = "wavelets",
                          nlevels = -1),
               "'nlevels' must be a single positive integer or NULL.")
  expect_error(fitFPLS_DA(data = spectra100, method = "wavelets",
                          nlevels = 1.5),
               "'nlevels' must be a single positive integer or NULL.")

  # level
  expect_error(fitFPLS_DA(data = spectra100, method = "wavelets", level = -1),
               "'level' must be a single positive integer.")
  expect_error(fitFPLS_DA(data = spectra100, method = "wavelets", level = 0),
               "'level' must be a single positive integer.")
  expect_error(fitFPLS_DA(data = spectra100, method = "wavelets", level = 1.5),
               "'level' must be a single positive integer.")
})
