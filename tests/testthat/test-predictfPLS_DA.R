# Load built-in data required for testing
data("spectra100", package = "maldiscrim")

test_that("predict.FPLS_DA works correctly for fitted values and new data", {
  # Fit a quick model to use for predictions
  model <- fitFPLS_DA(data = spectra100, method = "bsplines", nbasis = 50, nfolds = 2)

  # Test prediction when omitting 'newdata' (Fitted values / Training data)
  pred_fitted <- predict(model)

  expect_type(pred_fitted, "list")
  expect_named(pred_fitted, c("class", "probability", "lda_res"))
  expect_equal(length(pred_fitted$class), nrow(spectra100))
  expect_equal(length(pred_fitted$probability), nrow(spectra100))

  # 3. Test prediction with explicit 'newdata'
  new_spectra <- spectra100[1:5, ]
  pred_new <- predict(model, newdata = new_spectra)

  expect_vector(pred_new$class, "factor")
  expect_equal(length(pred_new$class), 5)
})

test_that("predict.FPLS_DA handles edge cases and raises expected errors", {
  model <- fitFPLS_DA(data = spectra100, method = "bsplines", nbasis = 50, nfolds = 2)

  # Check that providing an empty or improperly formatted dataframe raises an error
  empty_data <- data.frame()
  expect_error(predict(model, newdata = empty_data))
})
