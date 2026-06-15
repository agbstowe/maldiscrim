#' Fit a Functional Neural Network (FNN) for MALDI-TOF Spectrum Classification
#'
#' @description
#' `fitFNN` trains a Functional Neural Network on preprocessed MALDI-TOF mass spectra for strain classification.
#' The model is built from two functional convolutional layers followed by a functional dense output layer,
#' and is fitted via the Python TensorFlow/Keras backend through \code{reticulate}.
#'
#' @param data A numeric matrix Rows are biological replicates; row names must be the strain labels.
#' @param trainSize Numeric. Proportion of \code{data} used for training; the remaining proportion is held out internally
#' as a validation set for early stopping (monitored via \code{monitor}). Default is \code{0.8}.
#' @param batchSize Integer. Number of samples per gradient update. Default is \code{32}.
#' @param nbEpochs Integer. Maximum number of training epochs. Default is \code{50}.
#' @param stepsPerEpoch Integer. Number of steps per epoch. Default is \code{16}.
#' @param patience Integer. Number of epochs with no improvement before early stopping. Default is \code{5}.
#' @param monitor Character. Metric monitored for early stopping; one of
#'   \code{"loss"}, \code{"val_loss"}, \code{"accuracy"}, \code{"val_accuracy"}. Default is \code{"val_loss"}.
#' @param firstLayer A named list specifying the first convolutional layer with fields:
#'   \code{nFilters} (integer, default \code{50}),
#'   \code{nFunctions} (integer, default \code{100}),
#'   \code{resolution} (integer, default \code{500}),
#'   \code{basisType} (character, either \code{"Fourier"} or \code{"Legendre"}, default \code{"Fourier"}),
#'   \code{activation} (character, default \code{"elu"}).
#' @param secondLayer A named list specifying the second convolutional layer with fields:
#'   \code{nFilters} (integer, default \code{5}),
#'   \code{nFunctions} (integer, default \code{70}),
#'   \code{resolution} (integer, default \code{200}),
#'   \code{basisType} (character, either \code{"Fourier"} or \code{"Legendre"}, default \code{"Fourier"}),
#'   \code{activation} (character, default \code{"elu"}).
#' @param finalLayer A named list specifying the final dense layer with fields:
#'   \code{nFunctions} (integer, default \code{50}),
#'   \code{resolution} (integer, default \code{19966}),
#'   \code{basisType} (character, either \code{"Fourier"} or \code{"Legendre"}, default \code{"Fourier"}),
#'   \code{activation} (character, default \code{"softmax"}),
#'   \code{pooling} (logical, default \code{TRUE}).
#' @param outputPath Character. Path to the directory where the trained model and training outputs will be saved.
#' @param verbose Logical. If \code{TRUE}, prints training progress and messages. Default is \code{TRUE}.
#'
#' @details
#' \code{fitFNN} internally splits \code{data} according to \code{trainSize} to create a validation set used for early stopping.
#' The full model returned is the one obtained after this internal training procedure; \code{predict.FNN}
#' is then used to evaluate the model on new, independent data.
#'
#' @return It returns an object of class \code{"FNN"}. This object is a named list containing:
#' \item{modelPath}{Path to the saved Keras model file.}
#' \item{history}{A named list of training metrics (loss, accuracy) per epoch.}
#' \item{classNames}{Character vector of strain names in classification order.}
#' \item{nClass}{Number of distinct strain classes.}
#' \item{nSample}{Total number of samples used.}
#' \item{trainSize}{Proportion of data used for training.}
#' \item{batchSize}{Batch size used during training.}
#' \item{nbEpochs}{Maximum number of epochs set.}
#' \item{patience}{Early stopping patience used.}
#' \item{monitor}{Metric monitored for early stopping.}
#' \item{firstLayer}{List of first convolutional layer parameters.}
#' \item{secondLayer}{List of second convolutional layer parameters.}
#' \item{finalLayer}{List of final dense layer parameters.}
#' \item{outputPath}{Path where model outputs were saved.}
#'
#' @seealso \code{\link{ProcessMALDI}}, \code{\link{predict.FNN}}, \code{\link{summary.FNN}}
#'
#' @examples
#' \donttest{
#' # data is the output of ProcessMALDI()
#'
#' # Fit with default parameters
#' fnn_model <- fitFNN(data = spectra100)
#'
#' # Fit with custom architecture
#' fnn_model <- fitFNN(
#'   data        = spectra100,
#'   nbEpochs    = 100,
#'   patience    = 10,
#'   firstLayer  = list(nFilters = 64, nFunctions = 120,
#'                      resolution = 500, basisType = "Fourier", activation = "elu"),
#'   secondLayer = list(nFilters = 8,  nFunctions = 80,
#'                      resolution = 200, basisType = "Fourier", activation = "elu"),
#'   finalLayer  = list(nFunctions = 50, resolution = 19966,
#'                      basisType = "Fourier", activation = "softmax", pooling = TRUE)
#' )
#' }
#'
#' @export
fitFNN <- function(data, trainSize  = 0.8, batchSize  = 32L, nbEpochs  = 50L, stepsPerEpoch = 16L,
                   patience     = 5L, monitor      = "val_loss",
                   firstLayer   = list(nFilters = 50L, nFunctions = 100L,
                                       resolution = 500L, basisType = "Fourier", activation = "elu"),
                   secondLayer  = list(nFilters = 5L,  nFunctions = 70L,
                                       resolution = 200L, basisType = "Fourier", activation = "elu"),
                   finalLayer   = list(nFunctions = 50L, resolution = 19966L,
                                       basisType = "Fourier", activation = "softmax", pooling = TRUE),
                   outputPath   = "outputs/fnn", verbose  = TRUE) {

  # Input validation
  if (!is.matrix(data) && !is.data.frame(data)) {
    stop("'data' must be a numeric matrix or data.frame (output of ProcessMALDI).")
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

  # Layer parameter validation
  .validateLayerOptions(firstLayer,  c("nFilters", "nFunctions", "resolution", "basisType", "activation"), "firstLayer")
  .validateLayerOptions(secondLayer, c("nFilters", "nFunctions", "resolution", "basisType", "activation"), "secondLayer")
  .validateLayerOptions(finalLayer,  c("nFunctions", "resolution", "basisType", "activation", "pooling"),  "finalLayer")

  # Python environment check
  if (!reticulate::py_available()) {
    stop("Python is not available. Please run maldiscrim_install_python() first.")
  }

  # Export data to CSV
  nClass     <- length(levels(as.factor(rownames(data))))
  csvPath    <- file.path(outputPath, "FNN_input.csv")

  meta <- .exportForFNN(data = data, filePath = csvPath, verbose = verbose)

  if (verbose) message("Starting FNN training...")

  # Python training call
  result <- .callPythonFNN(
    csvPath       = meta$filePath,
    nClass        = meta$nClass,
    trainSize     = trainSize,
    batchSize     = batchSize,
    nbEpochs      = nbEpochs,
    stepsPerEpoch = stepsPerEpoch,
    patience      = patience,
    monitor       = monitor,
    firstLayer    = firstLayer,
    secondLayer   = secondLayer,
    finalLayer    = finalLayer,
    outputPath    = outputPath,
    verbose       = verbose
  )

  if (verbose) message("FNN training complete.")

  # Output object
  modelObj <- list(
    modelPath   = result$modelPath,
    history     = result$history,
    classNames  = meta$classNames,
    nClass      = meta$nClass,
    nSample     = meta$nSample,
    trainSize   = trainSize,
    batchSize   = batchSize,
    nbEpochs    = nbEpochs,
    patience    = patience,
    monitor     = monitor,
    firstLayer  = firstLayer,
    secondLayer = secondLayer,
    finalLayer  = finalLayer,
    outputPath  = outputPath
  )

  class(modelObj) <- "FNN"
  return(modelObj)
}
