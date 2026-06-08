# Load package built-in data required for the baseline simulation input
data("spectra100", package = "maldiscrim")

test_that("SimulateSpectrum outputs correct dimensions and names", {
  # Run the simulation for n = 15 synthetic spectra
  target_n <- 15
  sim_result <- SimulateSpectrum(data = spectra100, n = target_n, k = 3, plot = FALSE)

  # Check output structure
  expect_true(is.matrix(sim_result))
  expect_equal(nrow(sim_result), target_n)
  expect_equal(ncol(sim_result), ncol(spectra100))

  # Check naming conventions based on cluster strains and original features
  expect_equal(colnames(sim_result), colnames(spectra100))
  expect_match(rownames(sim_result)[1], "strain [1-3]")
})

test_that("SimulateSpectrum graphical workflow operates without failure", {
  pdf(file = tempfile())
  on.exit(dev.off())

  expect_no_error(SimulateSpectrum(data = spectra100, n = 5, k = 2, plot = TRUE))
})

test_that("SimulateSpectrum handles edge cases and parameters boundaries", {
  # Check k default fallback behavior when NULL is supplied
  sim_k_null <- SimulateSpectrum(data = spectra100, n = 6, k = NULL, plot = FALSE)
  expect_equal(nrow(sim_k_null), 6)

  # Check invalid input safeguards (empty dataset)
  empty_matrix <- matrix(numeric(0), nrow = 0, ncol = 0)
  expect_error(SimulateSpectrum(data = empty_matrix, n = 10))
})
