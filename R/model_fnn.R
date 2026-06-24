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
#' \item{dataPath}{Path to the CSV file used for training.}
#' \item{history}{A named list of training metrics (loss, accuracy) per epoch.}
#' \item{nParams}{Integer. Total number of trainable parameters in the network.}
#' \item{labels}{The original factor vector of strain labels assigned to each row of \code{data}.}
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
  if (!is.numeric(batchSize) || length(batchSize) != 1 ||
      batchSize < 1 || batchSize != round(batchSize))
    stop("'batchSize' must be a single positive integer.")

  if (!is.numeric(stepsPerEpoch) || length(stepsPerEpoch) != 1 ||
      stepsPerEpoch < 1 || stepsPerEpoch != round(stepsPerEpoch))
    stop("'stepsPerEpoch' must be a single positive integer.")

  # convert value as integer for python use
  batchSize     <- as.integer(batchSize)
  nbEpochs      <- as.integer(nbEpochs)
  stepsPerEpoch <- as.integer(stepsPerEpoch)
  patience      <- as.integer(patience)

  firstLayer$nFilters    <- as.integer(firstLayer$nFilters)
  firstLayer$nFunctions  <- as.integer(firstLayer$nFunctions)
  firstLayer$resolution  <- as.integer(firstLayer$resolution)
  secondLayer$nFilters   <- as.integer(secondLayer$nFilters)
  secondLayer$nFunctions <- as.integer(secondLayer$nFunctions)
  secondLayer$resolution <- as.integer(secondLayer$resolution)
  finalLayer$nFunctions  <- as.integer(finalLayer$nFunctions)
  finalLayer$resolution  <- as.integer(finalLayer$resolution)

  # Layer parameter validation
  .validateLayerOptions(firstLayer,  c("nFilters", "nFunctions", "resolution", "basisType", "activation"), "firstLayer")
  .validateLayerOptions(secondLayer, c("nFilters", "nFunctions", "resolution", "basisType", "activation"), "secondLayer")
  .validateLayerOptions(finalLayer,  c("nFunctions", "resolution", "basisType", "activation", "pooling"),  "finalLayer")

  # Python environment check
  # if (!reticulate::py_available()) {
  #   stop("Python is not available. Please run maldiscrim_install_python() first.")
  # }

  # Python environment check
  .checkPythonEnv()

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
    dataPath    = meta$filePath,
    history     = result$history,
    nParams     = result$nParams,
    labels      = as.factor(rownames(data)),
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


#' Predict Method for FNN Model Fits
#'
#' @description
#' Obtains predictions (classified strains and posterior probabilities) from a fitted \code{FNN} model object.
#'
#' @method predict FNN
#'
#' @param object A fitted model object of class \code{"FNN"}.
#' @param newdata An optional matrix or data frame of new, MALDI-TOF spectral intensities to predict.
#' If omitted (\code{NULL}), the function retrieves predictions on the data used to fit the model.
#' @param threshold Numeric. A classification probability threshold (between 0 and 1).
#' If the maximum softmax probability for a sample is below this threshold,
#' its predicted class is labeled as \code{"Doubty"}. Default is \code{NULL}.
#' @param ... Further arguments passed to or from other methods.
#'
#' @details
#' When new spectral data is supplied via \code{newdata}, it is assumed to be unlabeled and is
#' exported to a temporary CSV file containing only spectral intensities, then passed to the trained Keras model.
#' The softmax output layer of the network provides, for each observation, a probability for every strain class;
#' the predicted class is the one with the highest probability.
#'
#' @return A named list containing three components:
#' \item{class}{A character vector of predicted strain assignments for each observation.}
#' \item{probability}{A numeric vector indicating the maximum softmax probability associated with each prediction.}
#' \item{softmax}{A numeric matrix of softmax probabilities for every class, one row per observation.}
#'
#' @seealso \code{\link{fitFNN}}, \code{\link{summary.FNN}}
#'
#' @examples
#' \donttest{
#' # data is the output of ProcessMALDI()
#' fnn_model <- fitFNN(data = spectra100)
#'
#' # Predict on the training data by omitting 'newdata'
#' predict_res <- predict(fnn_model)
#' table(predict_res$class)
#'
#' # Example with 'newdata' on a subset of spectra100
#' new_spectra <- spectra100[1:5, ]
#' predict_res <- predict(fnn_model, newdata = new_spectra)
#' print(predict_res$class)
#' print(predict_res$probability)
#' }
#'
#' @importFrom stats predict
#' @export
predict.FNN <- function(object, newdata = NULL, threshold = NULL, ...) {

  if (!inherits(object, "FNN")) {
    stop("'object' must be a fitted FNN model.")
  }
  if (!is.null(threshold) && (!is.numeric(threshold) || threshold < 0 || threshold > 1)) {
    stop("'threshold' must be a single numeric value between 0 and 1.")
  }

  object <- .resolveFNNPaths(object)

  # Python environment check
  .checkPythonEnv()

  if (is.null(newdata)) {
    predictPath <- object$dataPath
    hasLabels   <- TRUE

  } else {
    if (!is.matrix(newdata) && !is.data.frame(newdata)) {
      stop("'newdata' must be a matrix or a data.frame.")
    }
    tmpCsv      <- tempfile(pattern = "fnn_newdata_", fileext = ".csv")
    .writeFeaturesCSV(data = newdata, filePath = tmpCsv)
    predictPath <- tmpCsv
    hasLabels   <- FALSE
  }

  # Prediction
  predResult <- .callPythonPredict(
    modelPath   = object$modelPath,
    dataPath    = predictPath,
    nClass      = object$nClass,
    batchSize   = object$batchSize,
    hasLabels   = hasLabels
  )

  softmaxProbs <- predResult$softmax
  classNames   <- object$classNames

  maxProbs    <- apply(softmaxProbs, 1, max)
  finalClass  <- classNames[apply(softmaxProbs, 1, which.max)]

  if (!is.null(threshold)) {
    finalClass[maxProbs < threshold] <- "Doubty"
  }

  return(list(
    class       = finalClass,
    probability = maxProbs,
    softmax     = softmaxProbs
  ))
}


#' Summary Method for FNN Model Fits
#'
#' @description
#' Prints a summary of a fitted \code{FNN} model, covering the network architecture,
#' training trajectory, and in-sample classification performance.
#'
#' @method summary FNN
#'
#' @param object A fitted model object of class \code{"FNN"}.
#' @param ... Further arguments passed to or from other methods.
#'
#' @details
#' The summary is organised into three sections:
#' \itemize{
#'   \item **Architecture** — filters, basis type, number of basis functions, resolution and activation for each of the
#'     two convolutional layers and the final dense layer, along with the total number of trainable parameters.
#'   \item **Training** — epochs run versus epochs set, the best epoch retained,
#'     early stopping settings, and a train/validation comparison of loss and accuracy at that best epoch.
#'   \item **Performance** — number of samples and groups, in-sample confusion matrix, overall accuracy, and per-group recall.
#' }
#'
#' In-sample predictions are obtained by calling \code{predict.FNN} on the training data (i.e. \code{newdata = NULL}).
#' For the full layer-by-layer Keras architecture, use \code{architectureFNN(object)} instead.
#'
#' @return Invisibly returns a named list with the following elements, allowing programmatic access to the computed summary statistics:
#' \item{nParams}{Integer. Total number of trainable parameters.}
#' \item{nbEpochsRun}{Integer. Number of epochs actually run before stopping.}
#' \item{bestEpoch}{Integer. Epoch index whose weights were retained by early stopping.}
#' \item{trainLoss}{Numeric. Training loss at the best epoch.}
#' \item{trainAccuracy}{Numeric. Training accuracy at the best epoch.}
#' \item{valLoss}{Numeric. Validation loss at the best epoch.}
#' \item{valAccuracy}{Numeric. Validation accuracy at the best epoch.}
#' \item{nSamples}{Integer. Total number of training samples.}
#' \item{nGroups}{Integer. Number of distinct groups.}
#' \item{confusion}{Table. In-sample confusion matrix.}
#' \item{accuracy}{Numeric. Overall in-sample accuracy (0-1).}
#' \item{recall}{Named numeric vector. Per-group recall (0-1).}
#'
#' @seealso \code{\link{fitFNN}}, \code{\link{predict.FNN}}, \code{\link{architectureFNN}}
#'
#' @examples
#' \donttest{
#' # data is the output of ProcessMALDI()
#' fnn_model <- fitFNN(data = spectra100)
#' summary(fnn_model)
#'
#' # Programmatic access to summary statistics
#' s <- summary(fnn_model)
#' s$accuracy
#' s$confusion
#' }
#'
#' @export
summary.FNN <- function(object, ...) {

  object <- .resolveFNNPaths(object)

  # In-sample predictions
  pred        <- predict.FNN(object, newdata = NULL)
  true_labels <- as.character(object$labels)
  # confusion   <- table(Predicted = pred$class, Actual = true_labels)
  cf          <- caret::confusionMatrix(as.factor(pred$class), as.factor(true_labels))
  confusion   <- cf$table
  accuracy    <- sum(diag(confusion)) / sum(confusion)

  # Per-group recall
  recall <- sapply(colnames(confusion), function(g) {
    confusion[g, g] / sum(confusion[, g])
  })

  # Best epoch: early stopping restores weights from the epoch with the best monitored value.
  monitorIsLoss <- grepl("loss", object$monitor, fixed = TRUE)
  monitorValues <- object$history[[object$monitor]]
  bestEpoch     <- if (monitorIsLoss) which.min(monitorValues) else which.max(monitorValues)
  nbEpochsRun   <- length(object$history$loss)

  trainLoss     <- object$history$loss[[bestEpoch]]
  trainAccuracy <- object$history$accuracy[[bestEpoch]]
  valLoss       <- object$history$val_loss[[bestEpoch]]
  valAccuracy   <- object$history$val_accuracy[[bestEpoch]]

  sep_thick <- strrep("=", 56)
  sep_thin  <- strrep("-", 56)

  cat(sep_thick, "\n")
  cat(" FNN Model Summary\n")
  cat(sep_thick, "\n\n")

  # 1: Architecture
  cat("[ 1 ] Architecture\n")
  cat(sep_thin, "\n")
  cat(sprintf("  Layer 1 (conv)  : %d filters, %s (%d fn, res %d), activation : %s\n",
              object$firstLayer$nFilters, object$firstLayer$basisType,
              object$firstLayer$nFunctions, object$firstLayer$resolution,
              object$firstLayer$activation))
  cat(sprintf("  Layer 2 (conv)  : %d filters, %s (%d fn, res %d), activation : %s\n",
              object$secondLayer$nFilters, object$secondLayer$basisType,
              object$secondLayer$nFunctions, object$secondLayer$resolution,
              object$secondLayer$activation))
  cat(sprintf("  Layer 3 (dense) : %d classes, %s (%d fn, res %d), activation : %s\n",
              object$nClass, object$finalLayer$basisType,
              object$finalLayer$nFunctions, object$finalLayer$resolution,
              object$finalLayer$activation))
  cat(sprintf("  Trainable parameters : %s\n", format(object$nParams, big.mark = ",", scientific = FALSE)))
  cat("\n")

  # 2: Training
  cat("\n[ 2 ] Training\n")
  cat(sep_thin, "\n")
  cat(sprintf("  Epochs run / set    : %d / %d   (best epoch: %d)\n", nbEpochsRun, object$nbEpochs, bestEpoch))
  cat(sprintf("  Early stopping      : monitor = %s, patience = %d\n", object$monitor, object$patience))
  cat(sprintf("  Batch size / split  : %d  |  %.0f%% train / %.0f%% validation\n",
              object$batchSize, object$trainSize * 100, (1 - object$trainSize) * 100))
  cat("\n")
  cat(sprintf("  %-14s%-12s%-12s\n", "", "Train", "Validation"))
  cat(sprintf("  %-14s%-12.4f%-12.4f\n", "Loss", trainLoss, valLoss))
  cat(sprintf("  %-14s%-12s%-12s\n", "Accuracy",
              sprintf("%.1f%%", trainAccuracy * 100), sprintf("%.1f%%", valAccuracy * 100)))
  cat("\n")

  # 3: Performance
  cat(sprintf("\n[ 3 ] Performance  (n = %d, %d groups)\n", length(true_labels), length(unique(true_labels))))
  cat(sep_thin, "\n")
  cat(sprintf("  Overall accuracy : %.1f%%\n\n", accuracy * 100))
  cat("  Confusion matrix :\n")
  print(confusion)
  cat("\n  Per-group recall:\n")
  recall_df <- data.frame(Group = names(recall), Recall = sprintf("%.1f%%", recall * 100), row.names = NULL)
  print(recall_df, row.names = FALSE, quote = FALSE)

  # Invisible return for programmatic access
  invisible(list(
    nParams       = object$nParams,
    nbEpochsRun   = nbEpochsRun,
    bestEpoch     = bestEpoch,
    trainLoss     = trainLoss,
    trainAccuracy = trainAccuracy,
    valLoss       = valLoss,
    valAccuracy   = valAccuracy,
    nSamples      = length(true_labels),
    nGroups       = length(unique(true_labels)),
    confusion     = confusion,
    accuracy      = accuracy,
    recall        = recall
  ))
}


#' Display the Keras Architecture of a Fitted FNN Model
#'
#' @description
#' Displays the full Keras architecture of a fitted \code{FNN} model, layer by layer.
#'
#' @param object A fitted model object of class \code{"FNN"}.
#'
#' @details
#' The trained Keras model is reloaded and its native architecture summary is printed, showing for each layer its name,
#' output shape, and number of parameters. For a concise statistical report on training and classification performance, use
#' \code{summary(object)} instead.
#'
#' @return Invisibly returns \code{object}.
#'
#' @seealso \code{\link{fitFNN}}, \code{\link{summary.FNN}}
#'
#' @examples
#' \donttest{
#' # data is the output of ProcessMALDI()
#' fnn_model <- fitFNN(data = spectra100)
#' architectureFNN(fnn_model)
#' }
#'
#' @export
architectureFNN <- function(object) {

  if (!inherits(object, "FNN")) {
    stop("'object' must be a fitted FNN model.")
  }

  object <- .resolveFNNPaths(object)

  .checkPythonEnv()

  fnn <- .loadFNNModel(object$modelPath)

  cat(strrep("=", 56), "\n")
  cat(" FNN Keras Architecture\n")
  cat(strrep("=", 56), "\n\n")

  fnn$summary(print_fn = function(line) cat(line, "\n"))

  invisible(object)
}
