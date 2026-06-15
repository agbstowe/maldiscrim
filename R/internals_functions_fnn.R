# These functions are internal to the FNN part of the package.
# They are not exported and not visible in the package documentation.

# ══════════════════════════════════════════════════════════════════════════════
# .exportForFNN
# ══════════════════════════════════════════════════════════════════════════════

#' Export preprocessed spectra to CSV for FNN ingestion
#'
#' @description
#' Converts the numeric matrix produced by \code{ProcessMALDI} into a CSV file formatted for direct ingestion by the Python FNN pipeline.
#' The CSV contains spectral intensities as predictors followed by one-hot encoded class labels derived from the row names of \code{data}.
#'
#' @param data A numeric matrix as returned by \code{ProcessMALDI}. Rows are biological replicates; row names must be the strain labels.
#' @param filePath Character. Full path (including file name and \code{.csv} extension) where the CSV will be written.
#' @param verbose Logical. If \code{TRUE}, prints a confirmation message with the output path and matrix dimensions.
#'
#' @return It returns a list with:
#' \item{filePath}{The path of the written CSV file.}
#' \item{nSample}{Total number of rows written.}
#' \item{nClass}{Number of distinct strain classes.}
#' \item{classNames}{Character vector of strain names in one-hot order.}
#'
#' @keywords internal
.exportForFNN <- function(data, filePath, verbose = TRUE) {

  # Input validation
  if (!is.matrix(data) && !is.data.frame(data)) {
    stop("'data' must be a numeric matrix or data.frame.")
  }

  if (is.null(rownames(data)) || any(rownames(data) == "")) {
    stop("'data' must have non-empty row names corresponding to strain labels.")
  }

  if (!is.character(filePath) || length(filePath) != 1 || nchar(filePath) == 0) {
    stop("'filePath' must be a single non-empty character string.")
  }

  # Create output directory
  out_dir <- dirname(filePath)
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
    if (verbose) message(sprintf("Directory created: %s", out_dir))
  }

  # One-hot encoding
  labels      <- as.factor(rownames(data))
  classNames <- levels(labels)
  nClass     <- length(classNames)
  nSample   <- nrow(data)

  # Build one-hot matrix
  onehot <- matrix(0L, nrow = nSample, ncol = nClass,
                   dimnames = list(NULL, classNames))

  for (i in seq_len(nSample)) {
    onehot[i, as.integer(labels[i])] <- 1L
  }

  # write CSV
  intensity_colnames        <- paste0("mz_", seq_len(ncol(data)))
  export_matrix             <- cbind(data, onehot)
  colnames(export_matrix)   <- c(intensity_colnames, classNames)

  utils::write.csv(export_matrix, file = filePath, row.names = FALSE)

  if (verbose) {
    message(sprintf(
      "FNN input CSV written to: %s\n  Samples : %d\n  Classes : %d ", filePath, nSample, nClass))
  }

  # return
  invisible(list(
    filePath    = filePath,
    nSample   = nSample,
    nClass     = nClass,
    classNames = classNames
  ))
}


# ══════════════════════════════════════════════════════════════════════════════
# .validateLayerOptions
# ══════════════════════════════════════════════════════════════════════════════

#' Validate layer option lists for FNN
#'
#' @description
#' Checks that a layer specification list contains all required fields with non-null values.
#'
#' @param layerList A named list of layer parameters.
#' @param requiredFields Character vector of required field names.
#' @param layerName Character. Name of the layer, used in error messages.
#'
#' @keywords internal
.validateLayerOptions <- function(layerList, requiredFields, layerName) {

  if (!is.list(layerList)) {
    stop(sprintf("'%s' must be a named list.", layerName))
  }
  missing <- setdiff(requiredFields, names(layerList))
  if (length(missing) > 0) {
    stop(sprintf("'%s' is missing required field(s): %s.",layerName, paste(missing, collapse = ", ")))
  }
  if ("basisType" %in% names(layerList) && !layerList$basisType %in% c("Fourier", "Legendre")) {
    stop(sprintf("'%s$basisType' must be either \"Fourier\" or \"Legendre\".", layerName))
  }
}


# ══════════════════════════════════════════════════════════════════════════════
# .callPythonFNN
# ══════════════════════════════════════════════════════════════════════════════

#' Call the Python FNN training pipeline via reticulate
#'
#' @description
#' Passes all training parameters to the Python FNN backend and returns the path to the saved model and the training history.
#'
#' @param csvPath Character. Path to the input CSV produced by \code{.exportForFNN}.
#' @param nClass Integer. Number of strain classes.
#' @param trainSize Numeric. Train/test split ratio.
#' @param batchSize Integer. Batch size.
#' @param nbEpochs Integer. Maximum number of training epochs.
#' @param stepsPerEpoch Integer. Steps per epoch.
#' @param patience Integer. Early stopping patience.
#' @param monitor Character. Metric monitored for early stopping.
#' @param firstLayer Named list of first convolutional layer parameters.
#' @param secondLayer Named list of second convolutional layer parameters.
#' @param finalLayer Named list of final dense layer parameters.
#' @param outputPath Character. Directory where model outputs are saved.
#' @param verbose Logical. If \code{TRUE}, prints progress messages.
#'
#' @return A named list with:
#' \item{modelPath}{Path to the saved Keras model file.}
#' \item{history}{Named list of training metrics per epoch.}
#'
#' @keywords internal
.callPythonFNN <- function(csvPath, nClass, trainSize, batchSize, nbEpochs, stepsPerEpoch, patience, monitor,
                           firstLayer, secondLayer, finalLayer, outputPath, verbose) {

  # Load Python FNN module
  pyPath <- system.file("python", package = "maldiscrim")
  dml    <- reticulate::import_from_path("maldiscrim_fnn.data_model_loading", path = pyPath)

  load_data   <- dml$load_data
  setup_model <- dml$setup_model

  tf          <- reticulate::import("tensorflow")
  keras       <- reticulate::import("keras")

  # Output folder
  trainedPath  <- file.path(outputPath, "trainedFNN")
  if (!dir.exists(trainedPath)) dir.create(trainedPath, recursive = TRUE)

  # Data loading
  datasets     <- load_data(csvPath, ratio = trainSize, n_labels = nClass, batch_size = batchSize)
  trainData    <- datasets[[1]]
  testData     <- datasets[[2]]

  inputShape   <- reticulate::py_eval(
    sprintf("(%d, 1)", reticulate::py_to_r(testData$element_spec[[1]]$shape[[2]]))
  )

  # Build filter and layer options
  filterOptions <- list(
    list(n_filters     = firstLayer$nFilters,
         basis_options = list(n_functions = firstLayer$nFunctions, resolution  = firstLayer$resolution, basis_type  = firstLayer$basisType),
         activation    = firstLayer$activation),
    list(n_filters     = secondLayer$nFilters,
         basis_options = list(n_functions = secondLayer$nFunctions, resolution  = secondLayer$resolution, basis_type  = secondLayer$basisType),
         activation    = secondLayer$activation)
  )

  layerOptions  <- list(
    list(n_neurons     = nClass,
         basis_options = list(n_functions = finalLayer$nFunctions, resolution  = finalLayer$resolution, basis_type  = finalLayer$basisType),
         activation    = finalLayer$activation, pooling = finalLayer$pooling)
  )

  # Model setup and training
  fnn      <- setup_model(inputShape, filterOptions, layerOptions)
  stopper  <- keras$callbacks$EarlyStopping(monitor = monitor, patience = as.integer(patience), restore_best_weights = TRUE)

  training <- fnn$fit(trainData,
    epochs           = as.integer(nbEpochs),
    steps_per_epoch  = as.integer(stepsPerEpoch),
    validation_data  = testData,
    callbacks        = list(stopper),
    verbose          = as.integer(verbose)
  )

  # Save model
  modelPath <- file.path(trainedPath, "fnn.keras")
  fnn$save(modelPath)
  if (verbose) message(sprintf("Model saved to: %s", modelPath))

  # Save training history
  history      <- reticulate::py_to_r(training$history)
  historyPath  <- file.path(trainedPath, "training_history.json")
  jsonlite::write_json(history, historyPath)
  if (verbose) message(sprintf("Training history saved to: %s", historyPath))

  # return
  invisible(list(
    modelPath = modelPath,
    history   = history
  ))
}
