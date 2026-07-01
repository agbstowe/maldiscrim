#' Transfer a Pre-Trained Functional Neural Network (FNN) to New MALDI-TOF Data
#'
#' @description
#' \code{transferFNN} performs transfer learning from a previously trained \code{"FNN"} model object onto a new, preprocessed
#' MALDI-TOF dataset. The new dataset may correspond to new strains of the same species used to train the source model,
#' or to a different species altogether; the number of target classes is free to differ from the source model.
#'
#' @param sourceModel A fitted model object of class \code{"FNN"}, used as the pre-trained source for transfer learning. If \code{NULL} (default),
#' the reference model bundled with the package is used automatically.
#' @param data A numeric matrix of new, target MALDI-TOF spectra. Rows are biological replicates; row names must be the (new) strain labels. May
#' contain new strains of the same species as \code{sourceModel}, strains of a different species, or both.
#' @param strategy Character. One approach of \code{"feature_extraction"} (default), \code{"fine_tune_partial"}, \code{"fine_tune_full"}.
#' Ignored if \code{nFreezeLayers} is supplied. See Details.
#' @param nFreezeLayers Integer or \code{NULL}. If supplied, overrides \code{strategy} and directly sets the number of leading functional layers
#' of \code{sourceModel} to freeze. Default is \code{NULL}.
#' @param finalLayer A named list specifying the new head (\code{FunctionalDense}) layer, with fields:
#'   \code{nFunctions} (integer, default \code{50}),
#'   \code{resolution} (integer; must match the spatial resolution of the source model's backbone output -- see \code{architectureFNN(sourceModel)}),
#'   \code{basisType} (character, either \code{"Fourier"} or \code{"Legendre"}, default \code{"Fourier"}),
#'   \code{activation} (character, default \code{"softmax"}),
#'   \code{pooling} (logical, default \code{TRUE}).
#' @param outputPath Character. Path to the directory where the new model and training outputs will be saved.
#' @inheritParams fitFNN
#'
#' @details
#' The source model is never modified: \code{transferFNN} loads it, truncates it just before its final \code{FunctionalDense} layer,
#' freezes a number of leading layers determined by \code{strategy} (or overridden via \code{nFreezeLayers}), and grafts a freshly initialized
#' \code{FunctionalDense} head configured via \code{finalLayer}, exactly as in \code{fitFNN} onto the truncated backbone. The resulting model is
#' then fitted on \code{data} using the same early-stopping logic as \code{fitFNN}.
#'
#' Three built-in strategies are available, expressed as the number of leading functional layers (both \code{FunctionalConvolution} layers, in
#' network order) of the source model that are frozen during training:
#' \itemize{
#'   \item \code{"feature_extraction"} (default): both \code{FunctionalConvolution} layers are frozen; only the new head is trained.
#'   Recommended when the new dataset is small, or when the new strains belong to the same species as the source model.
#'   \item \code{"fine_tune_partial"}: only the first \code{FunctionalConvolution} layer is frozen; the second one and the new head are trained.
#'   Recommended when transferring to a different species with a sufficiently large new dataset.
#'   \item \code{"fine_tune_full"}: no layer is frozen; the entire backbone and the new head are trained.
#' }
#' \code{nFreezeLayers} can be supplied to override the number of frozen layers directly, bypassing \code{strategy}.
#'
#' Importantly, \code{transferFNN} never embeds or redistributes the data used to train the source model:
#' the saved model object only contains learned weights, architecture, and class names, never the original
#' spectra. Only \code{data} (the new, target dataset supplied) is used to train the new head.
#'
#'
#' @return An object of class \code{"FNN"}, identical in structure to the output of \code{\link{fitFNN}}, with two additional fields for traceability:
#' \item{sourceModelPath}{Path to the source model's saved Keras model file.}
#' \item{strategy}{Character. The transfer strategy used.}
#' \item{nFreezeLayers}{Integer. The actual number of leading layers frozen.}
#'
#' @seealso \code{\link{fitFNN}}, \code{\link{predict.FNN}}, \code{\link{summary.FNN}}, \code{\link{architectureFNN}}
#'
#' @examples
#' \donttest{
#' # Default usage: the package's bundled pre-trained model is used as source.
#' transfer_model <- transferFNN(data = new_spectra, strategy = "feature_extraction")
#'
#' predict(transfer_model)
#' summary(transfer_model)
#'
#' # Advanced usage: supply a custom source model trained with fitFNN().
#' source_model  <- fitFNN(data = spectra100)
#' transfer_model <- transferFNN(sourceModel = source_model, data = new_spectra, strategy = "fine_tune_partial")
#' }
#'
#' @export
transferFNN <- function(data, sourceModel = NULL, strategy = c("feature_extraction", "fine_tune_partial", "fine_tune_full"),
                        nFreezeLayers = NULL, trainSize = 0.8, batchSize = 32L, nbEpochs = 50L, stepsPerEpoch = 16L,
                        patience = 5L, monitor = "val_loss",
                        finalLayer = list(nFunctions = 50L, resolution = 19966L, basisType = "Fourier", activation = "softmax", pooling = TRUE),
                        outputPath = "outputs/transferFNN", verbose = TRUE) {

  strategy <- match.arg(strategy)

  # Load bundled pre-trained model when sourceModel is not supplied
  if (is.null(sourceModel)) {
    bundledPath <- system.file("models", "fnn_base_model.h5", package = "maldiscrim")
    if (!nzchar(bundledPath)) {
      stop("No bundled pre-trained model found in inst/models/fnn_base_model.h5. ",
           "Please supply a fitted FNN object via the 'sourceModel' argument.")
    }
    sourceModel <- list(
      modelPath   = bundledPath,
      firstLayer  = NULL,
      secondLayer = NULL
    )
    class(sourceModel) <- "FNN"
    if (verbose) cli::cli_inform("Using the package's bundled pre-trained FNN model as source.")
  }

  # Input validation
  if (!inherits(sourceModel, "FNN")) {
    stop("'sourceModel' must be a fitted FNN model or NULL.")
  }
  if (!is.matrix(data) && !is.data.frame(data)) {
    stop("'data' must be a numeric matrix or data.frame.")
  }
  if (is.null(rownames(data)) || any(rownames(data) == "")) {
    stop("'data' must have non-empty row names corresponding to strain labels.")
  }
  if (!is.numeric(trainSize) || length(trainSize) != 1 || trainSize <= 0 || trainSize >= 1) {
    stop("'trainSize' must be a single numeric value in (0, 1).")
  }
  if (!is.numeric(nbEpochs) || length(nbEpochs) != 1 || nbEpochs < 1 || nbEpochs != round(nbEpochs)) {
    stop("'nbEpochs' must be a single positive integer.")
  }
  if (!is.numeric(patience) || length(patience) != 1 || patience < 1 || patience != round(patience)) {
    stop("'patience' must be a single positive integer.")
  }
  if (!is.character(monitor) || length(monitor) != 1 || !monitor %in% c("loss", "val_loss", "accuracy", "val_accuracy")) {
    stop('\'monitor\' must be one of "loss", "val_loss", "accuracy", "val_accuracy".')
  }
  if (!is.character(outputPath) || length(outputPath) != 1 || nchar(outputPath) == 0) {
    stop("'outputPath' must be a single non-empty character string.")
  }
  if (!is.numeric(batchSize) || length(batchSize) != 1 || batchSize < 1 || batchSize != round(batchSize))
    stop("'batchSize' must be a single positive integer.")
  if (!is.numeric(stepsPerEpoch) || length(stepsPerEpoch) != 1 || stepsPerEpoch < 1 || stepsPerEpoch != round(stepsPerEpoch))
    stop("'stepsPerEpoch' must be a single positive integer.")
  if (!is.null(nFreezeLayers) && (!is.numeric(nFreezeLayers) || length(nFreezeLayers) != 1 ||
                                  nFreezeLayers < 0 || nFreezeLayers != round(nFreezeLayers))) {
    stop("'nFreezeLayers' must be NULL or a single non-negative integer.")
  }

  # convert value as integer for python use
  batchSize     <- as.integer(batchSize)
  nbEpochs      <- as.integer(nbEpochs)
  stepsPerEpoch <- as.integer(stepsPerEpoch)
  patience      <- as.integer(patience)

  finalLayer$nFunctions <- as.integer(finalLayer$nFunctions)
  finalLayer$resolution <- as.integer(finalLayer$resolution)

  # Layer parameter validation (reuses the same validator as fitFNN)
  .validateLayerOptions(finalLayer, c("nFunctions", "resolution", "basisType", "activation", "pooling"), "finalLayer")

  # Python environment check
  .checkPythonEnv()

  # Resolve source model path
  sourceModel <- .resolveFNNPaths(sourceModel)

  # Resolve freeze strategy
  resolvedNFreezeLayers <- .nFreezeLayersForStrategy(strategy, nFreezeLayers)
  resolvedStrategy      <- if (is.null(nFreezeLayers)) strategy else "custom"

  # Export new data to CSV
  nClass  <- length(levels(as.factor(rownames(data))))
  csvPath <- file.path(outputPath, "transferFNN_input.csv")

  meta <- .exportForFNN(data = data, filePath = csvPath, verbose = verbose)

  if (verbose) {
    cli::cli_alert_warning("This process may take longer because an automatic prediction step will follow the training.")
    cli::cli_inform(sprintf("Starting FNN transfer learning (strategy: %s, frozen layers: %d)...", resolvedStrategy, resolvedNFreezeLayers))
  }

  # Python transfer learning call
  result <- .callPythonTransferFNN(
    sourceModelPath = sourceModel$modelPath,
    csvPath         = meta$filePath,
    nClass          = meta$nClass,
    nFreezeLayers   = resolvedNFreezeLayers,
    trainSize       = trainSize,
    batchSize       = batchSize,
    nbEpochs        = nbEpochs,
    stepsPerEpoch   = stepsPerEpoch,
    patience        = patience,
    monitor         = monitor,
    finalLayer      = finalLayer,
    outputPath      = outputPath,
    verbose         = verbose
  )

  if (verbose) cli::cli_alert_success("FNN transfer learning complete!")

  # Output object
  modelObj <- list(
    modelPath       = result$modelPath,
    dataPath        = meta$filePath,
    history         = result$history,
    nParams         = result$nParams,
    labels          = as.factor(rownames(data)),
    classNames      = meta$classNames,
    nClass          = meta$nClass,
    nSample         = meta$nSample,
    trainSize       = trainSize,
    batchSize       = batchSize,
    nbEpochs        = nbEpochs,
    patience        = patience,
    monitor         = monitor,
    firstLayer      = sourceModel$firstLayer,
    secondLayer     = sourceModel$secondLayer,
    finalLayer      = finalLayer,
    outputPath      = outputPath,
    sourceModelPath = sourceModel$modelPath,
    strategy        = resolvedStrategy,
    nFreezeLayers   = resolvedNFreezeLayers
  )

  class(modelObj) <- "FNN"

  modelObj$fittedValues <- predict(modelObj)

  return(modelObj)
}
