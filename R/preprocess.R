#' Import and Preprocessing of Biological Replicates from MALDI-TOF Strains
#'
#' @description
#' `ProcessMALDI` automates the complete preprocessing pipeline for raw MALDI-TOF mass spectra organized in strain-specific subdirectories.
#' It handles variance stabilization, smoothing, baseline removal, intensity calibration, technical
#' replicate alignment, and computes Mean Mass Spectra (MSP) for each biological replicate.
#'
#' @param dataFile Character. Path to the parent directory containing strain subfolders with raw spectra files.
#' @param methodBaseline Character. The baseline subtraction method to be passed to `MALDIquant::removeBaseline`. Default is `"SNIP"`.
#' @param halfWindowSize Integer. The half-window size used for Savitzky-Golay intensity smoothing. Default is `10`.
#'
#' @details
#' The function scans the `dataFile` directory for subfolders, where each subfolder
#' represents a unique biological strain. For each strain, the following sequential pipeline is executed via the `MALDIquant` ecosystem:
#' \enumerate{
#'   \item **Importation:** Reads raw data files using \code{\link[MALDIquantForeign]{import}}.
#'   \item **Variance Stabilization:** Applies a square root (\code{sqrt}) transformation.
#'   \item **Smoothing:** Smooths intensities using the Savitzky-Golay algorithm.
#'   \item **Baseline Correction:** Subtracts the background baseline using the SNIP algorithm (with 100 iterations).
#'   \item **Normalization:** Calibrates intensities using Total Ion Current (\code{TIC}).
#'   \item **Grouping & Alignment:** Groups spectra into biological replicates by combining
#'   the metadata sample full name and the acquisition month. Technical replicates within each group
#'   are aligned using a Lowess warping method.
#'   \item **Average Spectrum:** Computes a Mean Mass Spectrum (MSP) for each biological replicate.
#' }
#' To ensure a consistent matrix structure for downstream functional data analysis (FDA) and deep learning (CNN),
#' all intensity vectors are standardized to a fixed length of 20,664 data points.
#' Missing values (\code{NA}) introduced during padding are replaced with 0.
#'
#' @return A numeric matrix where each row represents a processed biological replicate
#' and columns represent the fixed-length spectral intensities (20,664 variables at all).
#'
#' @seealso \code{\link[MALDIquantForeign]{import}}, \code{\link[MALDIquant]{alignSpectra}},
#' \code{\link[MALDIquant]{averageMassSpectra}}
#'
#' @examples
#' \dontrun{
#' # Assuming you have a directory structure like:
#' # "data/Strain_A/...", "data/Strain_B/..." , ...
#'
#' path <- "path/to/data"
#'
#' # Process all raw spectra into a standardized intensity matrix
#' spectra_matrix <- ProcessMALDI(dataFile = path, methodBaseline = "SNIP", halfWindowSize = 10)
#'
#' # Preview result
#' dim(spectra_matrix)
#' head(rownames(spectra_matrix))
#' }
#'
#' # Note: The pre-processed dataset included in this package (e.g., 'spectra100')
#' # was built using this exact pipeline.
#'
#'
#' @export
ProcessMALDI <- function(dataFile, methodBaseline = "SNIP", halfWindowSize = 10) { # length intensity adding

  strains <- list.dirs(dataFile, full.names = FALSE, recursive = FALSE)
  rep_bio <- do.call(rbind, lapply(strains, function(strain_){
    file_strain <- file.path(dataFile, strain_)

    ## Importation and Preprocessing
    spectra <- MALDIquantForeign::import(file_strain, verbose = FALSE)
    # Variance Stabilization
    spectra <- MALDIquant::transformIntensity(spectra, method = "sqrt")
    # Smoothing
    spectra <- MALDIquant::smoothIntensity(spectra, method = "SavitzkyGolay", halfWindowSize = halfWindowSize)
    # Baseline removing
    spectra <- MALDIquant::removeBaseline(spectra, method = methodBaseline, iterations = 100)
    # Normalization
    spectra <- MALDIquant::calibrateIntensity(spectra, method = "TIC")

    # Grouping by biologic replicas with replica names and month
    strain <- factor(sapply(spectra, function(x) {
      paste0(MALDIquant::metaData(x)$fullName, "_", lubridate::month(MALDIquant::metaData(x)$acquisitionDate))
    }))
    strain_spectra <- split(spectra, strain)

    # Alignment and MSP computation
    msp_list <- lapply(strain_spectra, function(spectralist) {
      aligned <- MALDIquant::alignSpectra(spectralist, halfWindowSize = 20, SNR = 2, tolerance = 0.001, warpingMethod = "lowess")
      return(MALDIquant::averageMassSpectra(aligned, method = "mean"))
    })

    # msp_strain <- do.call(cbind, lapply(msp_list, function(msp_obj) {
    #   intensity <- msp_obj@intensity
    #   length(intensity) <- 20664
    #   intensity[is.na(intensity)] <- 0
    #   return(intensity)
    # }))

    msp_strain <- vapply(msp_list, function(msp_obj) {
      intensity <- msp_obj@intensity
      length(intensity) <- 20664
      intensity[is.na(intensity)] <- 0
      return(intensity)
    }, FUN.VALUE = numeric(20664))

    msp_strain <- t(msp_strain)

    rownames(msp_strain) <- rep(strain, nrow(msp_strain))

    return(msp_strain)
  }))

  return(rep_bio)
}
