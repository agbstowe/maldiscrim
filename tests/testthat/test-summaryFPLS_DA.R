# Load precomputed model for all summary tests
model_bsplines <- fitFPLS_DA(precomputed = TRUE)


test_that("summary.FPLS_DA returns a list with correct names", {
  s <- summary(model_bsplines)
  expect_type(s, "list")
  expect_named(s, c("method", "ncompOpt", "nRemoved", "varExp", "nSamples", "nGroups", "Frequency", "confusion", "accuracy", "recall"))
})

test_that("summary.FPLS_DA fields have correct types and valid ranges", {
  s <- summary(model_bsplines)

  # Types
  expect_type(s$method,   "character")
  expect_type(s$ncompOpt, "integer")
  expect_type(s$nRemoved, "integer")
  expect_type(s$varExp,   "double")
  expect_type(s$nSamples, "integer")
  expect_type(s$nGroups,  "integer")
  expect_s3_class(s$Frequency, "table")
  expect_s3_class(s$confusion, "table")
  expect_type(s$accuracy, "double")
  expect_type(s$recall,   "double")

  # Ranges
  expect_true(s$ncompOpt >= 1)
  expect_true(s$nRemoved >= 0)
  expect_true(s$varExp   >  0 && s$varExp <= 100)
  expect_true(s$nSamples >  0)
  expect_true(s$nGroups  >= 2)
  expect_true(s$accuracy >= 0 && s$accuracy <= 1)
  expect_true(all(s$recall >= 0 & s$recall <= 1))
})


test_that("summary.FPLS_DA statistics are internally consistent", {
  s <- summary(model_bsplines)

  expect_equal(s$nSamples, sum(s$confusion))
  expect_equal(sum(s$Frequency), s$nSamples)
  expect_equal(s$nGroups, ncol(s$confusion))
  expect_equal(s$nGroups, length(s$recall))
  expected_acc <- sum(diag(s$confusion)) / sum(s$confusion)
  expect_equal(s$accuracy, expected_acc)
  expect_equal(s$ncompOpt, model_bsplines$ncompOpt)
  expect_equal(s$nRemoved, sum(model_bsplines$const_idx))
})


test_that("summary() dispatches correctly on FPLS_DA objects", {
  expect_no_error(summary(model_bsplines))
  s <- summary(model_bsplines)
  expect_false(is.null(s))
})


test_that("summary.FPLS_DA prints expected sections to console", {
  output <- capture.output(summary(model_bsplines))
  combined <- paste(output, collapse = "\n")

  expect_match(combined, "FPLS-DA Model Summary")
  expect_match(combined, "Functional Decomposition")
  expect_match(combined, "PLS Dimension Reduction")
  expect_match(combined, "Training Data")
  expect_match(combined, "In-Sample Classification Performance")
  expect_match(combined, "Overall accuracy")
  expect_match(combined, "Confusion matrix")
  expect_match(combined, "Per-group recall")
})
