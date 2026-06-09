#' Fit Functionnal Partiel Least Square and Discriminant Analysis (fPLS_DA) Model
#'
#' @description
#' `fitFPLS_DA` performs a supervised classification of MALDI-TOF mass spectra by combining functional data analysis (FDA)
#' dimensionality reduction techniques with Partial Least Squares Regression (PLSR) and Linear Discriminant Analysis (LDA).
#'
#'
#' @param data A numeric matrix of processed Mean Mass Spectra (MSP) where rows
#' represent biological replicates (with strain names as rownames) and columns represent spectral intensities.
#' @param method Character. The mathematical basis used for spectral decomposition: either `"bsplines"`
#' (B-spline basis functions) or `"wavelets"` (Discrete Wavelet Transform).
#' @param nbasis Integer. Number of basis functions used when `method = "bsplines"`. Default is `1050`.
#' @param ncomp Integer. Number of PLS components to retain. If `NULL` (default), the optimal number of
#' components is automatically determined using the one-sigma heuristic based on Cross-Validation RMSEP.
#' @param argvals A numeric vector specifying the evaluation grid (m/z indices) for the B-spline projection.
#' @param rangeval A numeric vector of length 2 defining the domain range for the B-spline basis.
#' @param nfolds Integer. Number of folds for the internal cross-validation performed during PLS fitting. Default is `5`.
#' @param filter Character. The wavelet filter type passed to \code{\link[wavelets]{dwt}} when `method = "wavelets"`. Default is `"la14"`.
#' @param boundary Character. The boundary handling method passed to \code{\link[wavelets]{dwt}}. Default is `"periodic"`.
#' @param nlevels Integer. Total number of decomposition levels for the wavelet transform.
#' If `NULL`, it is automatically set to the maximum possible level.
#' @param level Integer. The specific wavelet decomposition resolution level from which
#' coefficients are extracted for downstream classification. Default is `3`.
#' @param precomputed Logical. If `TRUE`, the function bypasses all heavy computations and returns a built-in
#' pre-trained model object. Primarily used to speed up examples and unit testing. Default is `FALSE`.
#'
#'
#' @details
#' The pipeline executes three sequential phases to achieve classification of functional spectral data:
#' \enumerate{
#'   \item **Functional Decomposition:** Depending on `method`, the high-dimensional spectral matrix is projected into either:
#'     \itemize{
#'       \item *B-Splines:* Continuous functional representations are built using a B-spline basis matrix via `fda::Data2fd`.
#'       \item *Wavelets:* Multi-resolution analysis is applied via the Discrete Wavelet Transform (`wavelets::dwt`). A built-in guardrail automatically adjusts `level` if the remaining coefficients are insufficient for modeling.
#'     }
#'   \item **PLS Dimension Reduction:** A Kernel Partial Least Squares Regression (`pls::plsr`) is fitted using dummy-coded groups as responses.
#'   \item **Classification:** Low-variance components are filtered out, and a Linear Discriminant Analysis (`MASS::lda`) is fitted on the remaining PLS scores to generate the final decision boundaries.
#' }
#'
#' @return An object of class `"FPLS_DA"`. This object is a named list containing:
#' \item{PLSmodel}{The fitted PLS model object returned by `pls::plsr`.}
#' \item{LDAmodel}{The fitted LDA model object returned by `MASS::lda`.}
#' \item{method}{A character string indicating the functional decomposition method used.}
#' \item{basis}{The B-spline basis object (if `method = "bsplines"`), otherwise `NULL`.}
#' \item{wavelets}{A list containing raw DWT outputs for each spectrum (if `method = "wavelets"`), otherwise `NULL`.}
#' \item{ncompOpt}{The optimised or user-defined number of PLS components retained.}
#' \item{const_idx}{A logical vector indicating which PLS components were excluded due to near-zero variance.}
#' \item{levels}{The structural names/categories of the groups (strains).}
#' \item{labels}{The original factor vector of group memberships assigned to each row.}
#' \item{argvals}{The evaluation grid used for B-spline projection, or `NULL`.}
#' \item{filter}{The wavelet filter name used, or `NULL`.}
#' \item{nlevels}{The effective wavelet decomposition level count, or `NULL`.}
#' \item{level}{The specific wavelet resolution level extracted, or `NULL`.}
#' \item{nbasis}{The total number of B-spline basis functions used, or `NULL`.}
#' \item{boundary}{The wavelet boundary method applied, or `NULL`.}
#'
#' @seealso \code{\link[pls]{plsr}}, \code{\link[MASS]{lda}}, \code{\link[fda]{Data2fd}}, \code{\link[wavelets]{dwt}}
#'
#' @examples
#' # Model loading using the precomputed cache
#' model_res <- fitFPLS_DA(precomputed = TRUE)
#' # print(model_res)
#'
#' \donttest{
#' # Real execution example using the package built-in data 'spectra100'
#' data("spectra100", package = "maldiscrim")
#'
#' # Fit with B-Splines functional basis
#' bspline_model <- fitFPLS_DA(
#'   data = spectra100, method = "bsplines", nbasis = 1050, nfolds = 5
#'   )
#'
#' # Fit with Wavelets multi-resolution decomposition
#' # wavelet_model <- fitFPLS_DA(
#' # data = spectra100, method = "wavelets", level = 3, nfolds = 5)
#' }
#'
#' @importFrom stats sd aggregate
#' @export
fitFPLS_DA <- function(data, method = c("bsplines", "wavelets"), nbasis = 1050, ncomp = NULL, argvals = 2001:22664,
                          rangeval = c(1999,22664), nfolds = 5,
                          filter = "la14", boundary = "periodic", nlevels = NULL, level = 3, precomputed = FALSE) {

  if (precomputed) {
    message("Loading precomputed FPLS_DA model from cache")
    if (exists("fpls_model", envir = asNamespace("maldiscrim"))) {
      return(get("fpls_model", envir = asNamespace("maldiscrim")))
    } else {
      stop("Cached model 'fpls_model' not found in the package data.")
    }
  }

  method <- match.arg(method)
  groups <- as.factor(rownames(data))

  # Input validation
  if (!is.matrix(data) && !is.data.frame(data)) {
    stop("'data' must be a matrix or a data.frame.")
  }
  if (!is.numeric(nbasis) || length(nbasis) != 1 || nbasis < 1 || nbasis != round(nbasis)) {
    stop("'nbasis' must be a single positive integer.")
  }
  if (!is.null(ncomp) && (!is.numeric(ncomp) || length(ncomp) != 1 || ncomp < 1 || ncomp != round(ncomp))) {
    stop("'ncomp' must be a single positive integer or NULL.")
  }
  if (!is.numeric(nfolds) || length(nfolds) != 1 || nfolds < 2 || nfolds != round(nfolds)) {
    stop("'nfolds' must be a single integer >= 2.")
  }
  if (!is.null(nlevels) && (!is.numeric(nlevels) || length(nlevels) != 1 || nlevels < 1 || nlevels != round(nlevels))) {
    stop("'nlevels' must be a single positive integer or NULL.")
  }
  if (!is.numeric(level) || length(level) != 1 || level < 1 || level != round(level)) {
    stop("'level' must be a single positive integer.")
  }

  # Basis decomposition
  res_fdaDecomposition <- .fdaDecomposition(data = data, method = method, nbasis = nbasis, rangeval = rangeval,
                              argvals  = argvals, filter  = filter, boundary = boundary, nlevels  = nlevels, level  = level)
  coef_matrix          <- res_fdaDecomposition$coef_matrix

  # PLS model
  y_dummy    <- fastDummies::dummy_cols(groups, remove_selected_columns = TRUE)
  max_search <- min(nrow(data) - 1L, 20L)
  pls_mod    <- pls::plsr(as.matrix(y_dummy) ~ coef_matrix, ncomp = max_search, method = "kernelpls",
                       validation = "CV", segments = nfolds)

  # Optimal number of components
  if (!is.null(ncomp)) {
    ncomp   <- min(ncomp, max_search)
    message(sprintf("User-supplied ncomp bounded to %d (maximum allowed for this dataset).", ncomp))
  } else {
    ncomp   <- .selectOptimalNcomp(pls_mod)
    message(sprintf("Optimal number of PLS components automatically set to %d.", ncomp))
  }

  # Score extraction and near-zero-variance filtering
  scores        <- pls_mod$scores[, seq_len(ncomp), drop = FALSE]
  group_means   <- aggregate(scores, list(groups), mean)[, -1]
  f1            <- apply(scores - as.matrix(group_means[groups, ]), 2, stats::sd)
  const_idx     <- f1 < 0.0001
  scores_clean  <- scores[, !const_idx, drop = FALSE]

  # LDA model
  lda_mod     <- MASS::lda(x = scores_clean, grouping = groups)

  # Output object
  model_obj <- list(
    pls_model = pls_mod,
    lda_model = lda_mod,
    method = method,
    basis = if(method == "bsplines") res_fdaDecomposition$basis_obj else NULL,
    wavelets = if(method == "wavelets") res_fdaDecomposition$wt_output else NULL,
    ncompOpt = ncomp,
    const_idx = const_idx,
    levels = levels(groups),
    labels = groups,
    argvals = if(method == "bsplines") argvals else NULL,
    filter = if(method == "wavelets") filter else NULL,
    nlevels = if(method == "wavelets") res_fdaDecomposition$nlevels else NULL,
    level = if(method == "wavelets") res_fdaDecomposition$level else NULL,
    nbasis = if(method == "bsplines") nbasis else NULL,
    boundary = if(method == "wavelets") boundary else NULL
  )

  class(model_obj) <- "FPLS_DA"
  return(model_obj)
}


#' Predict Method for FPLS_DA Model Fits
#'
#' @description
#' Obtains predictions (classified strains and posterior probabilities) from a fitted \code{FPLS_DA} model object.
#'
#' @method predict FPLS_DA
#'
#' @param object A fitted model object of class \code{"FPLS_DA"}.
#' @param newdata An optional matrix or data frame of new MALDI-TOF spectral intensities to predict.
#' If omitted (\code{NULL}), the function automatically retrieves the coefficients from the training data and returns the fitted predictions.
#' @param threshold Numeric. A classification probability threshold (between 0 and 1).
#' If the maximum posterior probability for a sample is below this threshold, its predicted class is labeled as \code{"Doubty"}.
#' Default is \code{NULL} (no threshold applied).
#' @param ... Further arguments passed to or from other methods.
#' @importFrom stats predict
#'
#' @details
#' When new spectral data is supplied via \code{newdata}, the function projects it into the exact same functional domain
#' space established during training (using either the original B-spline basis parameters or the specific wavelet filter and resolution level).
#'
#' The extracted functional coefficients are then projected onto the optimal Partial Least Squares (PLS) components.
#' Finally, after removing near-zero variance components, the Linear Discriminant Analysis boundary profiles
#' evaluate the scores to output the most probable strain class along with its associated posterior probability.
#'
#' @return A named list containing three components:
#' \item{class}{A character vector of predicted class/strain assignments for each observation.}
#' \item{probability}{A numeric vector indicating the maximum posterior probability (confidence score) associated with each prediction.}
#' \item{lda_res}{The raw list output from the underlying \code{\link[MASS]{predict.lda}} step.}
#'
#' @seealso \code{\link{fitFPLS_DA}}, \code{\link[MASS]{predict.lda}}, \code{\link[pls]{predict.mvr}}
#'
#'
#' @examples
#' # Load the pre-trained example model from the package cache
#' model_res <- fitFPLS_DA(precomputed = TRUE)
#'
#' # Predict on the training data by omitting 'newdata' (Fitted values)
#' predict_res <- predict(model_res)
#' table(predict_res$class)
#'
#' \donttest{
#' # Example with 'newdata' using the built-in dataset
#' data("spectra100", package = "maldiscrim")
#'
#' # Let's simulate new data by taking a subset of spectra100
#' new_spectra <- spectra100[1:5, ]
#'
#' # Predict
#' predict_res <- predict(model_res, newdata = new_spectra )
#' print(predict_res$class)
#' print(predict_res$probability)
#' }
#'
#' @importFrom stats predict
#' @importFrom MASS lda
#' @import pls
#' @export
predict.FPLS_DA <- function(object, newdata = NULL, threshold = NULL, ...) {

  if (is.null(newdata)) {
    new_coefs <- object$pls_model$model[[2]]

  } else {
    if (!is.matrix(newdata) && !is.data.frame(newdata)) {
      stop("newdata must be a matrix or a data.frame.")
    }
    new_coefs <- .fdaDecomposeNew(newdata, object)
  }

  # PLS
  pls_scores <- predict(object$pls_model, comps = seq_len(object$ncompOpt), newdata = as.matrix(new_coefs), type = "scores")
  # Near-zero-variance filtering
  clean_scores           <- data.frame(pls_scores[, !object$const_idx])
  colnames(clean_scores) <- colnames(object$lda_model$means)
  # LDA
  lda_res     <- predict(object$lda_model, newdata = clean_scores)
  final_class <- as.character(lda_res$class)
  max_probs   <- apply(lda_res$posterior, 1, max)

  if (!is.null(threshold)) {
    final_class[max_probs < threshold] <- "Doubty"
  }

  return(list(
    class = final_class,
    probability = max_probs,
    lda_res = lda_res
  ))
}




#' Summary Method for FPLS_DA Model Fits
#'
#' @description
#' Prints a structured summary of a fitted \code{FPLS_DA} model, covering the decomposition settings, PLS dimension reduction,
#' training data composition, and in-sample classification performance.
#'
#' @method summary FPLS_DA
#'
#' @param object A fitted model object of class \code{"FPLS_DA"}.
#' @param ... Further arguments passed to or from other methods (currently ignored).
#'
#' @details
#' The summary is organised into three sections:
#' \itemize{
#'   \item **Model configuration** — functional decomposition method and its associated parameters (number of B-spline basis functions,
#'      or wavelet filter, boundary, decomposition level and total levels).
#'   \item **PLS dimension reduction** — number of components searched, retained, and removed by the near-zero variance filter, along with the
#'     cumulative variance explained by the retained components.
#'   \item **Training data and performance** — number of samples, group  distribution table, in-sample confusion matrix, overall accuracy, and
#'     per-group recall.
#' }
#'
#' In-sample predictions are obtained by calling \code{predict.FPLS_DA} on the training coefficients (i.e. \code{newdata = NULL}).
#'
#' @return Invisibly returns a named list with the following elements, allowing programmatic access to the computed summary statistics:
#' \item{method}{Character. Functional decomposition method.}
#' \item{ncompOpt}{Integer. Number of PLS components retained.}
#' \item{nRemoved}{Integer. Number of components removed by variance filter.}
#' \item{varExp}{Numeric. Cumulative variance explained by retained components (percentage).}
#' \item{nSamples}{Integer. Total number of training samples.}
#' \item{nGroups}{Integer. Number of distinct groups.}
#' \item{Frequency}{Table. Sample counts per group.}
#' \item{confusion}{Table. In-sample confusion matrix.}
#' \item{accuracy}{Numeric. Overall in-sample accuracy (0–1).}
#' \item{recall}{Named numeric vector. Per-group recall (0–1).}
#'
#' @seealso \code{\link{fitFPLS_DA}}, \code{\link{predict.FPLS_DA}}
#'
#' @examples
#' model <- fitFPLS_DA(precomputed = TRUE)
#' summary(model)
#'
#' # Programmatic access to summary statistics
#' s <- summary(model)
#' s$accuracy
#' s$confusion
#'
#' @export
summary.FPLS_DA <- function(object, ...) {

  # In-sample predictions
  pred        <- predict.FPLS_DA(object, newdata = NULL)
  true_labels <- as.character(object$labels)
  confusion   <- table(Predicted = pred$class, Actual = true_labels)
  accuracy    <- sum(diag(confusion)) / sum(confusion)

  # Per-group recall:
  recall <- sapply(colnames(confusion), function(g) {
    confusion[g, g] / sum(confusion[, g])
  })

  # Variance explained by retained components
  pct      <- object$pls_model$Xvar / object$pls_model$Xtotvar * 100
  var_exp  <- sum(pct[seq_len(object$ncompOpt)])
  nRemoved <- sum(object$const_idx)

  # Group distribution
  Frequency <- table(Freq = true_labels)

  sep_thick <- strrep("=", 55)
  sep_thin  <- strrep("-", 55)

  cat(sep_thick, "\n")
  cat(" FPLS-DA Model Summary\n")
  cat(sep_thick, "\n\n")

  # 1: Model configuration
  cat("[ 1 ] Functional Decomposition\n")
  cat(sep_thin, "\n")
  cat(sprintf("  Method           : %s\n", object$method))

  if (object$method == "bsplines") {
    cat(sprintf("  Number of bases  : %d\n", object$nbasis))
    cat(sprintf("  Argument range   : [%g, %g]\n", object$argvals[1], utils::tail(object$argvals, 1)))
  } else {
    cat(sprintf("  Wavelet filter   : %s\n",   object$filter))
    cat(sprintf("  Boundary         : %s\n",   object$boundary))
    cat(sprintf("  Resolution level : %d\n",   object$level))
    cat(sprintf("  Total levels     : %d\n",   object$nlevels))
  }

  # 2: PLS dimension reduction
  cat("\n[ 2 ] PLS Dimension Reduction\n")
  cat(sep_thin, "\n")
  cat(sprintf("  Components searched  : %d\n", object$pls_model$ncomp))
  cat(sprintf("  Components retained  : %d\n", object$ncompOpt))
  cat(sprintf("  Components removed   : %d  (near-zero variance)\n", nRemoved))
  cat(sprintf("  Cumulative variance  : %.1f%%\n", var_exp))

  # 3: Training data and performance
  cat("\n[ 3 ] Training Data\n")
  cat(sep_thin, "\n")
  cat(sprintf("  Samples : %d\n", length(true_labels)))
  cat(sprintf("  Groups  : %d\n\n", length(Frequency)))
  cat("  Distribution per group:\n")
  print(Frequency, quote = FALSE)

  cat("\n[ 4 ] In-Sample Classification Performance\n")
  cat(sep_thin, "\n")
  cat(sprintf("  Overall accuracy : %.1f%%\n\n", accuracy * 100))
  cat("  Confusion matrix :\n")
  print(confusion)
  cat("\n  Per-group recall:\n")
  recall_df <- data.frame(Group  = names(recall), Recall = sprintf("%.1f%%", recall * 100), row.names = NULL)
  print(recall_df, row.names = FALSE, quote = FALSE)
  cat("\n", sep_thick, "\n", sep = "")

  # Invisible return for programmatic access
  invisible(list(
    method      = object$method,
    ncompOpt    = object$ncompOpt,
    nRemoved   = nRemoved,
    varExp = var_exp,
    nSamples   = length(true_labels),
    nGroups    = length(Frequency),
    Frequency  = Frequency,
    confusion   = confusion,
    accuracy    = accuracy,
    recall      = recall
  ))
}
