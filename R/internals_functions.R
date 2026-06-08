# These functions are not exported and not visible in the package documentation.

# .fdaDecomposition

#' Functional decomposition of spectral data (internal)
#'
#' Projects a raw spectral matrix into a reduced coefficient matrix using either a B-spline basis or a Discrete Wavelet Transform (DWT).
#'
#' @param data A numeric matrix of spectral intensities (rows = samples, columns = m/z values).
#' @param method Character. Either `"bsplines"` or `"wavelets"`.
#' @param nbasis Integer. Number of B-spline basis functions. Used only when  `method = "bsplines"`.
#' @param rangeval Numeric vector of length 2. Domain range for the B-spline basis. Used only when `method = "bsplines"`.
#' @param argvals Numeric vector. Evaluation grid (m/z indices) for B-spline projection. Used only when `method = "bsplines"`.
#' @param filter Character. The wavelet filter type passed to \code{\link[wavelets]{dwt}} when `method = "wavelets"`.
#' @param boundary Character. The boundary handling method passed to \code{\link[wavelets]{dwt}}. Default is `"periodic"`.
#' @param nlevels Integer or `NULL`. Total number of wavelet decomposition levels. If `NULL`, set to the maximum possible level. Used only when `method = "wavelets"`.
#' @param level Integer. Wavelet resolution level from which coefficients are extracted. Used only when `method = "wavelets"`.
#'
#' @return A named list with the following elements:
#' \item{coef_matrix}{Numeric matrix of decomposition coefficients (rows = samples).}
#' \item{basis_obj}{The B-spline basis object ([fda::create.bspline.basis()]), or `NULL` if `method = "wavelets"`.}
#' \item{wt_output}{A list of raw [wavelets::dwt()] objects (one per sample), or `NULL` if `method = "bsplines"`.}
#' \item{nlevels}{The effective number of wavelet decomposition levels used or `NULL` if `method = "bsplines"`.}
#' \item{level}{The effective wavelet resolution level used, or `NULL` if `method = "bsplines"`.}
#'
#' @keywords internal
.fdaDecomposition <- function(data, method, nbasis, rangeval, argvals, filter, boundary, nlevels, level) {

  if (method == "bsplines") {
    message("Using B-spline functional basis decomposition.")
    basis_obj    <- fda::create.bspline.basis(rangeval = rangeval, nbasis = nbasis)
    fd_obj       <- fda::Data2fd(argvals   = argvals, y = t(data), basisobj = basis_obj)
    coef_matrix  <- t(fd_obj$coefs)

    return(list(
      coef_matrix = coef_matrix,
      basis_obj   = basis_obj,
      wt_output   = NULL,
      nlevels     = NULL,
      level       = NULL
    ))

  } else {
    message("Using Discrete Wavelet Transform (DWT) decomposition.")
    max_levels <- floor(log2(ncol(data)))

    # Resolve nlevels
    if (is.null(nlevels)) {
      nlevels <- max_levels
    } else if (nlevels > max_levels) {
      warning(sprintf("nlevels = %d exceeds the maximum allowed value. Reset to %d.", nlevels, max_levels))
      nlevels <- max_levels
    }
    # Theoretical number of coefficients at `level` ~ ncol(data) / 2^level
    n_coef_theoric <- floor(ncol(data) / (2^level))

    if (n_coef_theoric < 20) {
      new_level <- max(1L, floor(log2(ncol(data) / 20)))
      warning(sprintf(
        paste("level = %d yields too few wavelet coefficients for PLS", "(theoretical count: %d).
              Automatically adjusted to level = %d."), level, n_coef_theoric, new_level))
      level <- new_level
    }

    if (level > nlevels) {
      warning(sprintf(
        "level = %d exceeds nlevels = %d. Adjusted to nlevels.", level, nlevels))
      level <- nlevels
    }

    # Compute DWT for each sample and Extract wavelet detail coefficients at the chosen level
    wt_list <- apply(data, 1, function(x) {
      wavelets::dwt(x, filter = filter, n.levels = nlevels, boundary = boundary)
    })
    coef_matrix <- t(sapply(wt_list, function(wt) wt@W[[level]]))

    return(list(
      coef_matrix = coef_matrix,
      basis_obj   = NULL,
      wt_output   = wt_list,
      nlevels     = nlevels,
      level       = level
    ))
  }
}


# .selectOptimalNcomp

#' Select the optimal number of PLS components via one-sigma heuristic (internal)
#'
#' Determines the optimal number of PLS components by applying the one-sigma
#' rule to the cross-validated Root Mean Square Error of Prediction (RMSEP):
#' the smallest number of components whose RMSEP is within one standard error
#' of the minimum RMSEP is selected. The final value is the maximum across all response variables.
#'
#' @param pls_mod A fitted PLS model object returned by [pls::plsr()] with `validation = "CV"`.
#'
#' @return An integer giving the optimal number of PLS components.
#'
#' @keywords internal
.selectOptimalNcomp <- function(pls_mod) {

  rmsep_data  <- pls::RMSEP(pls_mod, estimate = "CV")$val
  n_responses <- dim(rmsep_data)[1]

  ncomp_vec <- sapply(seq_len(n_responses), function(i) {
    errors  <- drop(rmsep_data[i, 1, -1])
    min_idx <- which.min(errors)
    min_val <- errors[min_idx]
    # One-sigma threshold
    threshold <- min_val + stats::sd(errors) / sqrt(length(errors))
    which(errors <= threshold)[1L]
  })

  max(ncomp_vec)
}


# .fdaDecomposeNew

#' Project new spectral data into the functional space of a fitted model (internal)
#'
#' Reprojects new MALDI-TOF spectra into the exact same functional coefficient
#' space established during training, using either the stored B-spline basis or the original wavelet filter and resolution level.
#'
#' @param newdata A numeric matrix of new spectral intensities (rows = samples, columns = m/z values).
#' @param object A fitted model object of class `"FPLS_DA"`.
#'
#' @return A numeric matrix of functional coefficients (rows = samples) ready to be projected onto the PLS components.
#'
#' @keywords internal
.fdaDecomposeNew <- function(newdata, object) {

  if (object$method == "bsplines") {
    fd_new      <- fda::Data2fd(argvals  = object$argvals, y = t(newdata), basisobj = object$basis)
    coef_matrix <- t(fd_new$coefs)

  } else {
    coef_matrix <- t(apply(newdata, 1, function(x) {
      wt_x <- wavelets::dwt(x, filter  = object$filter, n.levels = object$nlevels, boundary = object$boundary)
      wt_x@W[[object$level]]
    }))
  }

  coef_matrix
}


# .computeLDAalpha

#' Compute LDA posterior probabilities and alpha transparency values (internal)
#'
#' Extracts the clean PLS scores from a fitted `"FPLS_DA"` object, runs [MASS::lda()] prediction, and maps each sample's maximum posterior
#' probability to an alpha transparency tier for plotting.
#'
#' @param object A fitted model object of class `"FPLS_DA"`.
#'
#' @return A numeric vector of alpha values, one per sample, following the three-tier rule: `1.0` if probability >= 0.8, `0.6` if >= 0.5, `0.2` otherwise.
#'
#' @keywords internal
.computeLDAalpha <- function(object) {

  raw_scores   <- object$pls_model$scores[, seq_len(object$ncompOpt),drop = FALSE]
  clean_scores <- data.frame(raw_scores[, !object$const_idx, drop = FALSE])
  colnames(clean_scores) <- colnames(object$lda_model$means)

  lda_res   <- stats::predict(object$lda_model, newdata = clean_scores)
  max_probs <- apply(lda_res$posterior, 1, max)

  ifelse(max_probs >= 0.8, 1.0,
         ifelse(max_probs >= 0.5, 0.6, 0.2))
}


# .getPalette

#' Select a colour palette adapted to the number of groups (internal)
#'
#' Returns a character vector of colours based on the number of strain groups and an optional user-supplied palette.
#'
#' @param n Integer. Number of distinct groups to colour.
#' @param palette Either `NULL` (automatic selection) or a user-supplied value accepted by [ggplot2::scale_color_manual()] / [plotly::plot_ly()] :
#'   a character vector of colours or a named `RColorBrewer` palette string.
#'
#' @return A character vector of `n` colours.
#'
#' @keywords internal
.getPalette <- function(n, palette = NULL) {

  if (!is.null(palette)) {
    return(palette)
  }

  if (n <= 8) {
    # RColorBrewer "Set1": 8 vivid, well-separated colours
    RColorBrewer::brewer.pal(n = max(3, n), name = "Set1")[seq_len(n)]
  } else {
    # pals::cols25(): 25 maximally distinct colours for large group counts
    if (n > 25) {
      stop(sprintf(
        "Number of groups (%d) exceeds the maximum supported (25). Please supply a custom palette via the 'palette' argument.",
        n
      ))
    }
    pals::cols25()[seq_len(n)]
  }
}
