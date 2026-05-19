#' Simulated MALDI-TOF Mass Spectrometry Dataset
#'
#' A simulated dataset containing 100 mass spectra generated using Functional Principal Components Analysis (FPCA).
#' This dataset serves as a standard reference benchmark for testing preprocessing, classification and clustering algorithms
#' within the package.
#'
#' @format A numeric matrix with 100 rows and 20,664 columns representing the mass-to-charge ratio (\eqn{m/z}) channels.
#' \describe{
#'   \item{Rows}{Simulated biological samples replicates.}
#'   \item{Columns}{Continuous \eqn{m/z} values where each cell contains the recorded intensity.}
#' }
#'
#' @details
#' The spectra were simulated by decomposing a source dataset of real raw mass spectra into a low-dimensional functional subspace.
#' Adaptive noise was injected into the estimated scores before reconstructing the signals onto the original dense
#' grid of 20,664 variables via a Karhunen-Loeve expansion.
#' This process preserves the continuous covariance structure and physical constraints typical of MALDI-TOF profiles.
#'
#' @usage data(spectra100)
#'
#' @source Simulated benchmark generated from real microbial MALDI-TOF source profiles.
#'
#' @examples
#' # Load the dataset
#' data(spectra100)
#'
#' # Check dimensions
#' dim(spectra100)
#'
#' # Visualise the first spectrum profile
#' plot(1:ncol(spectra100), spectra100[1, ], type = "l",
#'      xlab = "m/z", ylab = "Intensity")
"spectra100"


#' Precomputed fPLS_DA Model for Demonstration
#'
#' A pre-trained \code{fPLS_DA} model object fitted on the \code{spectra100} dataset
#' using B-splines decomposition with 1050 basis functions. This cached model
#' avoids redundant heavy computations during examples and vignettes execution.
#'
#' @format A list of class \code{fPLS_DA} containing PLS scores, loadings,
#' and the internal LDA model components.
#' @source Precomputed from the internal \code{spectra100} dataset.
"fpls_model"
