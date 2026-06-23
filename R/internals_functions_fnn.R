# These functions are internal to the FNN part of the package.
# They are not exported and not visible in the package documentation.

# ══════════════════════════════════════════════════════════════════════════════
# Python module cache
# TensorFlow and Keras register global classes at import time; re-importing
# them in the same Python session causes "already registered" errors.
# These helpers ensure each module is imported at most once per R session.
# ══════════════════════════════════════════════════════════════════════════════

.maldiscrimPyEnv <- new.env(parent = emptyenv())

.tf <- function() {
  if (is.null(.maldiscrimPyEnv$tf)) {
    .maldiscrimPyEnv$tf <- reticulate::import("tensorflow")
  }
  .maldiscrimPyEnv$tf
}

.keras <- function() {
  if (is.null(.maldiscrimPyEnv$keras)) {
    .maldiscrimPyEnv$keras <- reticulate::import("keras")
  }
  .maldiscrimPyEnv$keras
}

.dml <- function() {
  if (is.null(.maldiscrimPyEnv$dml)) {
    pyPath <- system.file("python", package = "maldiscrim")
    .maldiscrimPyEnv$dml <- reticulate::import_from_path(
      "maldiscrim_fnn.data_model_loading", path = pyPath
    )
  }
  .maldiscrimPyEnv$dml
}

.conv <- function() {
  if (is.null(.maldiscrimPyEnv$conv)) {
    pyPath <- system.file("python", package = "maldiscrim")
    .maldiscrimPyEnv$conv <- reticulate::import_from_path(
      "maldiscrim_fnn.convolution", path = pyPath
    )
  }
  .maldiscrimPyEnv$conv
}

.dense <- function() {
  if (is.null(.maldiscrimPyEnv$dense)) {
    pyPath <- system.file("python", package = "maldiscrim")
    .maldiscrimPyEnv$dense <- reticulate::import_from_path(
      "maldiscrim_fnn.dense", path = pyPath
    )
  }
  .maldiscrimPyEnv$dense
}


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
  # Requires inst/python/maldiscrim_fnn/ to be a proper Python package
  # (with an __init__.py) containing basis.py, convolution.py, dense.py
  # and data_model_loading.py, as these files use relative imports.
  dml    <- .dml()

  load_data   <- dml$load_data
  setup_model <- dml$setup_model

  tf    <- .tf()
  keras <- .keras()

  # Output folder
  trainedPath  <- file.path(outputPath, "trainedFNN")
  if (!dir.exists(trainedPath)) dir.create(trainedPath, recursive = TRUE)

  # Data loading
  datasets     <- load_data(csvPath, ratio = trainSize, n_labels = nClass, batch_size = batchSize)
  trainData    <- datasets[[1]]
  testData     <- datasets[[2]]

  # inputShape   <- reticulate::py_eval(
  #   sprintf("(%d, 1)", reticulate::py_to_r(testData$element_spec[[1]]$shape[[2]]))
  # )

  shape_list   <- reticulate::py_to_r(testData$element_spec[[1]]$shape$as_list())
  num_features <- shape_list[[2]]
  inputShape   <- as.integer(c(num_features, 1))

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

  # Total trainable parameters (computed once, while the model is in memory)
  nParams <- reticulate::py_to_r(fnn$count_params())

  # Save training history
  history      <- reticulate::py_to_r(training$history)
  historyPath  <- file.path(trainedPath, "training_history.json")
  jsonlite::write_json(history, historyPath)
  if (verbose) message(sprintf("Training history saved to: %s", historyPath))

  # return
  invisible(list(
    modelPath = modelPath,
    history   = history,
    nParams   = nParams
  ))
}


# ══════════════════════════════════════════════════════════════════════════════
# .writeFeaturesCSV
# ══════════════════════════════════════════════════════════════════════════════

#' Write spectral intensities to CSV without class labels
#'
#' @description
#' Writes a numeric matrix of spectral intensities to a CSV file, with no label columns.
#' Used for predicting on new, unlabeled spectra where the true class is unknown.
#'
#' @param data A numeric matrix of spectral intensities.
#' @param filePath Character. Full path where the CSV will be written.
#'
#' @keywords internal
.writeFeaturesCSV <- function(data, filePath) {

  out_dir <- dirname(filePath)
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }

  export_matrix           <- as.matrix(data)
  colnames(export_matrix) <- paste0("mz_", seq_len(ncol(export_matrix)))

  utils::write.csv(export_matrix, file = filePath, row.names = FALSE)
}


# ══════════════════════════════════════════════════════════════════════════════
# .loadFNNModel
# ══════════════════════════════════════════════════════════════════════════════

#' Load a trained Keras FNN model with its custom layers
#'
#' @description
#' Loads a saved \code{.keras} model file, registering the custom
#' \code{FunctionalConvolution} and \code{FunctionalDense} layers required to deserialize it.
#'
#' @param modelPath Character. Path to the saved \code{.keras} model file.
#'
#' @return The loaded Keras model object.
#'
#' @keywords internal
.loadFNNModel <- function(modelPath) {

  conv  <- .conv()
  dense <- .dense()
  keras <- .keras()

  customObjects <- list(
    FunctionalConvolution = conv$FunctionalConvolution,
    FunctionalDense        = dense$FunctionalDense
  )

  keras$models$load_model(modelPath, custom_objects = customObjects)
}


# ══════════════════════════════════════════════════════════════════════════════
# .loadFNNFeaturesTensor
# ══════════════════════════════════════════════════════════════════════════════

#' Load spectral features as a TensorFlow tensor
#'
#' @description
#' Reads the CSV file referenced by a fitted \code{FNN} object and returns
#' its spectral intensity columns, excluding the trailing one-hot label
#' columns, as a TensorFlow tensor ready for model inference.
#'
#' @details
#' The CSV is read with base R's \code{read.csv()}, which natively handles
#' both LF and CRLF line endings and the header row, then converted to a
#' TensorFlow tensor only at the end. This avoids low-level TensorFlow
#' string-slicing operations, whose Python-style negative-index syntax is
#' not reliably reproducible from plain R \code{:} expressions.
#'
#' @param object A fitted model object of class \code{"FNN"}.
#'
#' @return A TensorFlow tensor of shape \code{(n_samples, n_features, 1)},
#'   with the channel dimension already added.
#'
#' @keywords internal
.loadFNNFeaturesTensor <- function(object) {

  tf <- .tf()

  rawTable  <- utils::read.csv(object$dataPath, check.names = FALSE)
  nClass    <- length(object$classNames)
  nFeatures <- ncol(rawTable) - nClass

  featureMatrix <- as.matrix(rawTable[, seq_len(nFeatures), drop = FALSE])
  storage.mode(featureMatrix) <- "double"

  features <- tf$constant(featureMatrix, dtype = tf$float32)
  features <- tf$expand_dims(features, -1L)

  features
}


# ══════════════════════════════════════════════════════════════════════════════
# .callPythonPredict
# ══════════════════════════════════════════════════════════════════════════════

#' Call the Python FNN model for prediction via reticulate
#'
#' @description
#' Loads a trained Keras FNN model and produces softmax probabilities for the spectra contained in a given CSV file.
#'
#' @details
#' The CSV is read with base R's \code{read.csv()}, which natively handles
#' both LF and CRLF line endings and the header row, then converted to a
#' TensorFlow tensor only at the point of inference. This avoids low-level
#' TensorFlow string-slicing operations, whose Python-style negative-index
#' syntax is not reliably reproducible from plain R \code{:} expressions.
#'
#' @param modelPath Character. Path to the saved \code{.keras} model file.
#' @param dataPath Character. Path to the CSV file containing the spectra to predict.
#' @param nClass Integer. Number of strain classes.
#' @param batchSize Integer. Batch size used when feeding data to the model.
#' @param hasLabels Logical. If \code{TRUE}, the CSV contains \code{nClass} trailing one-hot label columns to be
#' dropped before prediction; if \code{FALSE}, the CSV contains spectral intensities only.
#'
#' @return A named list with:
#' \item{softmax}{A numeric matrix of softmax probabilities, one row per observation and one column per class.}
#'
#' @keywords internal
.callPythonPredict <- function(modelPath, dataPath, nClass, batchSize, hasLabels = TRUE) {

  tf <- .tf()
  fnn <- .loadFNNModel(modelPath)

  # Read and parse the CSV without train/test split.
  rawTable  <- utils::read.csv(dataPath, check.names = FALSE)
  nFeatures <- if (hasLabels) ncol(rawTable) - as.integer(nClass) else ncol(rawTable)

  featureMatrix <- as.matrix(rawTable[, seq_len(nFeatures), drop = FALSE])
  storage.mode(featureMatrix) <- "double"

  features <- tf$constant(featureMatrix, dtype = tf$float32)
  features <- tf$expand_dims(features, -1L)

  ds <- tf$data$Dataset$from_tensor_slices(features)
  ds <- ds$batch(as.integer(batchSize))

  # Getting Proba
  softmaxProbs <- reticulate::py_to_r(fnn$predict(ds))

  # Return
  invisible(list(
    softmax = softmaxProbs
  ))
}


# ══════════════════════════════════════════════════════════════════════════════
# .setMaldiscrimVenvHome
# ══════════════════════════════════════════════════════════════════════════════

#' Redirect the virtual environment storage location
#'
#' @description
#' Sets the \code{WORKON_HOME} environment variable to a package-specific
#' cache folder, so that the \code{"maldiscrim-env"} virtual environment is
#' always created and found in the same, safe location across sessions.
#'
#' @details
#' This avoids failures on machines where \code{reticulate}'s default virtual
#' environment folder sits inside a cloud-synced directory such as OneDrive.
#' Since \code{Sys.setenv()} only applies to the current R session, this
#' function must be called before any operation that creates, lists, or
#' activates the \code{"maldiscrim-env"} environment.
#'
#' @return Invisibly returns the path used as \code{WORKON_HOME}.
#'
#' @keywords internal
.setMaldiscrimVenvHome <- function() {

  venvHome <- tools::R_user_dir("maldiscrim", which = "cache")
  if (!dir.exists(venvHome)) dir.create(venvHome, recursive = TRUE)
  Sys.setenv(WORKON_HOME = venvHome)

  invisible(venvHome)
}


# ══════════════════════════════════════════════════════════════════════════════
# .checkPythonEnv
# ══════════════════════════════════════════════════════════════════════════════

#' Check that the maldiscrim Python environment is available and active
#'
#' @description
#' Verifies that the \code{"maldiscrim-env"} virtual environment exists, and
#' activates it for the current R session if it is not already the active
#' environment.
#'
#' @keywords internal
.checkPythonEnv <- function() {

  .setMaldiscrimVenvHome()

  envName <- "maldiscrim-env"

  if (!envName %in% reticulate::virtualenv_list()) {
    stop(sprintf(
      "The Python environment required for FNN is not set up.\nPlease run maldiscrim_install_python() once before using this function."
    ), call. = FALSE)
  }

  currentEnv <- tryCatch(reticulate::py_config()$pythonhome, error = function(e) "")

  if (!grepl(envName, currentEnv, fixed = TRUE)) {
    reticulate::use_virtualenv(envName, required = TRUE)
  }

  invisible(TRUE)
}


# ══════════════════════════════════════════════════════════════════════════════
# .extractNetworkSpace
# ══════════════════════════════════════════════════════════════════════════════

#' Extract Activation Spaces (Logits or Softmax) From Keras Model
#'
#' @keywords internal
.extractNetworkSpace <- function(object, on = "logits") {

  fnnModel <- .loadFNNModel(object$modelPath)
  tf       <- .tf()
  keras    <- .keras()

  # Read features using the same robust approach as .callPythonPredict
  rawTable  <- utils::read.csv(object$dataPath, check.names = FALSE)
  nClass    <- length(object$classNames)
  nFeatures <- ncol(rawTable) - nClass

  featureMatrix <- as.matrix(rawTable[, seq_len(nFeatures), drop = FALSE])
  storage.mode(featureMatrix) <- "double"

  features <- tf$constant(featureMatrix, dtype = tf$float32)
  features <- tf$expand_dims(features, -1L)

  if (on == "softmax") {
    return(reticulate::py_to_r(fnnModel$predict(features)))
  } else {
    # Reproduce the logits extraction from generate_results.py:
    # clone the model replacing the final softmax with a linear activation,
    # then copy the weights — this gives true pre-activation logits.
    clone_fn <- reticulate::py_func(function(layer) {
      config <- layer$get_config()
      if (!is.null(config[["activation"]]) && config[["activation"]] == "softmax") {
        config[["activation"]] <- "linear"
      }
      layer$`__class__`$from_config(config)
    })
    linearModel <- keras$models$clone_model(fnnModel, clone_function = clone_fn)
    linearModel$set_weights(fnnModel$get_weights())
    return(reticulate::py_to_r(linearModel$predict(features)))
  }
}


# ══════════════════════════════════════════════════════════════════════════════
# .loadFeaturesFromCSV
# ══════════════════════════════════════════════════════════════════════════════

#' Load and prepare feature tensor from the training CSV file
#'
#' Reads the CSV file stored at `object$dataPath`, separates the spectral
#' feature columns from the one-hot encoded class columns, and returns a
#' TensorFlow tensor ready for model inference.
#'
#' This helper centralises the data-loading logic shared by
#' `plotGradCamFNN`, `plotOcclusionFNN`, and `plotActivationFNN`, replacing
#' three identical copy-pasted blocks.
#'
#' @param object A fitted model object of class `"FNN"`.
#'
#' @return A named list with:
#' \item{features}{A TensorFlow float32 tensor of shape
#'   `[n_samples, n_features, 1]`, ready for direct model inference.}
#' \item{nClass}{Integer. Number of class columns in the CSV.}
#' \item{nFeatures}{Integer. Number of spectral feature columns.}
#'
#' @keywords internal
.loadFeaturesFromCSV <- function(object) {

  tf    <- .tf()
  nClass <- length(object$classNames)

  rawTable      <- utils::read.csv(object$dataPath, check.names = FALSE)
  nFeatures     <- ncol(rawTable) - nClass
  featureMatrix <- as.matrix(rawTable[, seq_len(nFeatures), drop = FALSE])
  storage.mode(featureMatrix) <- "double"

  features <- tf$constant(featureMatrix, dtype = tf$float32)
  features <- tf$expand_dims(features, -1L)

  list(
    features  = features,
    nClass    = nClass,
    nFeatures = nFeatures
  )
}
