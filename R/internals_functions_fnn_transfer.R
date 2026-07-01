# These functions are internal to the transfer learning part of the package.
# They are not exported and not visible in the package documentation.


# ══════════════════════════════════════════════════════════════════════════════
# .transfer
# ══════════════════════════════════════════════════════════════════════════════

# .transfer <- function() {
#   if (is.null(.maldiscrimPyEnv$transfer)) {
#     pyPath <- system.file("python", package = "maldiscrim")
#     .maldiscrimPyEnv$transfer <- reticulate::import_from_path(
#       "maldiscrim_fnn.transfer", path = pyPath
#     )
#   }
#   .maldiscrimPyEnv$transfer
# }


# ══════════════════════════════════════════════════════════════════════════════
# .nFreezeLayersForStrategy
# ══════════════════════════════════════════════════════════════════════════════

#' Resolve the number of layers to freeze for a given transfer strategy
#'
#' @description
#' Translates a named transfer learning \code{strategy} into the number of leading functional layers (in network order: both \code{FunctionalConvolution}
#' layers, then the original \code{FunctionalDense} head, which is always discarded and rebuilt) of the source model to freeze.
#'
#' @param strategy Character. One approach of \code{"feature_extraction"}, \code{"fine_tune_partial"}, \code{"fine_tune_full"}.
#' @param nFreezeLayers Integer or \code{NULL}. User-supplied override; if not \code{NULL}, takes precedence over \code{strategy}.
#'
#' @return Integer. Number of leading functional layers to freeze.
#'
#' @keywords internal
.nFreezeLayersForStrategy <- function(strategy, nFreezeLayers = NULL) {

  if (!is.null(nFreezeLayers)) {
    return(as.integer(nFreezeLayers))
  }

  switch(strategy,
    feature_extraction = 2L,
    fine_tune_partial   = 1L,
    fine_tune_full       = 0L,
    stop(sprintf("Unknown strategy: '%s'.", strategy))
  )
}


# ══════════════════════════════════════════════════════════════════════════════
# .callPythonTransferFNN
# ══════════════════════════════════════════════════════════════════════════════

#' Call the Python FNN transfer learning pipeline via reticulate
#'
#' @description
#' Loads a previously trained Keras FNN model, freezes its leading layers according to the chosen strategy, grafts a freshly initialized
#' \code{FunctionalDense} head onto the truncated backbone (allowing the number of target classes to differ from the source model), and fits the
#' resulting model on the new data.
#'
#' @param sourceModelPath Character. Path to the saved \code{.keras} model file of the source (pre-trained) model.
#' @param csvPath Character. Path to the input CSV produced by \code{.exportForFNN} for the new (target) data.
#' @param nClass Integer. Number of target strain classes.
#' @param nFreezeLayers Integer. Number of leading functional layers of the source model to freeze.
#' @param trainSize Numeric. Train/test split ratio.
#' @param batchSize Integer. Batch size.
#' @param nbEpochs Integer. Maximum number of training epochs.
#' @param stepsPerEpoch Integer. Steps per epoch.
#' @param patience Integer. Early stopping patience.
#' @param monitor Character. Metric monitored for early stopping.
#' @param finalLayer Named list of new head (\code{FunctionalDense}) layer parameters; same format as in \code{fitFNN}.
#' @param outputPath Character. Directory where model outputs are saved.
#' @param verbose Logical. If \code{TRUE}, prints progress messages.
#'
#' @return A named list with:
#' \item{modelPath}{Path to the saved Keras model file.}
#' \item{history}{Named list of training metrics per epoch.}
#' \item{nParams}{Integer. Total number of trainable parameters in the new model.}
#'
#' @keywords internal
.callPythonTransferFNN <- function(sourceModelPath, csvPath, nClass, nFreezeLayers, trainSize, batchSize, nbEpochs, stepsPerEpoch,
                                    patience, monitor, finalLayer, outputPath, verbose) {

  dml      <- .dml()
  transfer <- .transfer()

  load_data <- dml$load_data

  tf    <- .tf()
  keras <- .keras()

  # Output folder
  trainedPath <- file.path(outputPath, "trainedTransferFNN")
  if (!dir.exists(trainedPath)) dir.create(trainedPath, recursive = TRUE)

  # Load the source model
  baseModel <- .loadFNNModel(sourceModelPath)

  # Data loading (new/target data)
  datasets  <- load_data(csvPath, ratio = trainSize, n_labels = nClass, batch_size = batchSize)
  trainData <- datasets[[1]]
  testData  <- datasets[[2]]

  # New head (FunctionalDense) options -- same format as layer_options in setup_model
  newLayerOptions <- list(
    list(n_neurons     = nClass,
         basis_options = list(n_functions = finalLayer$nFunctions, resolution  = finalLayer$resolution, basis_type  = finalLayer$basisType),
         activation    = finalLayer$activation, pooling  = finalLayer$pooling)
  )

  # Build the transfer model: freeze source layers, graft new head
  fnn      <- transfer$setup_transfer_model(base_model = baseModel, n_freeze_layers = as.integer(nFreezeLayers), new_layer_options = newLayerOptions)
  stopper  <- keras$callbacks$EarlyStopping(monitor = monitor, patience = as.integer(patience), restore_best_weights = TRUE)

  training <- fnn$fit(trainData,
    epochs           = as.integer(nbEpochs),
    steps_per_epoch  = as.integer(stepsPerEpoch),
    validation_data  = testData,
    callbacks        = list(stopper),
    verbose          = as.integer(verbose)
  )

  # Save model
  modelPath <- file.path(trainedPath, "transfer_fnn.keras")
  fnn$save(modelPath)
  if (verbose) message(sprintf("Transfer model saved to: %s", modelPath))

  # Total trainable parameters (computed once, while the model is in memory)
  nParams <- reticulate::py_to_r(fnn$count_params())

  # Save training history
  history     <- reticulate::py_to_r(training$history)
  historyPath <- file.path(trainedPath, "transfer_training_history.json")
  jsonlite::write_json(history, historyPath)
  if (verbose) message(sprintf("Training history saved to: %s", historyPath))

  # return
  invisible(list(
    modelPath = modelPath,
    history   = history,
    nParams   = nParams
  ))
}
