# Load package built-in data required for the baseline simulation input
data("spectra100", package = "maldiscrim")

test_that("SimulateSpectrum outputs correct dimensions and names", {
  target_n   <- 15
  sim_result <- SimulateSpectrum(data = spectra100, n = target_n, k = 3,
                                 plot = FALSE)

  # Output structure
  expect_true(is.matrix(sim_result))
  expect_equal(nrow(sim_result), target_n)
  expect_equal(ncol(sim_result), ncol(spectra100))

  # Naming conventions
  expect_equal(colnames(sim_result), colnames(spectra100))
  # rownames follow "strain 1", "strain 2", "strain 3" pattern
  expect_match(rownames(sim_result)[1], "^strain [1-3]$")
})

test_that("SimulateSpectrum graphical workflow operates without failure", {
  pdf(file = tempfile())
  on.exit(dev.off())
  expect_no_error(
    SimulateSpectrum(data = spectra100, n = 5, k = 2, plot = TRUE)
  )
})

test_that("SimulateSpectrum handles edge cases and parameter boundaries", {
  # k = NULL should fall back to k = 3
  sim_k_null <- SimulateSpectrum(data = spectra100, n = 6, k = NULL,
                                 plot = FALSE)
  expect_equal(nrow(sim_k_null), 6)

  # Empty dataset should raise an error
  empty_matrix <- matrix(numeric(0), nrow = 0, ncol = 0)
  expect_error(SimulateSpectrum(data = empty_matrix, n = 10))
})
