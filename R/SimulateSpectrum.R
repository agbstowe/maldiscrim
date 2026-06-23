#' Simulate spectrum using Functionnal Principal Components Analysis (FPCA)
#'
#' @description
#' Simulates a set of n standard spectra based on real or simulated data using Functional Principal Components Analysis (FPCA) and adaptive noise injection.
#'
#' @param data A numeric matrix of Mass Spectrometry (MSP) data, where rows represent samples and columns represent m/z variables.
#' @param n Integer. Total number of spectra to simulate in the output matrix.
#' @param k Integer. Number of clusters/strains to identify in the latent space. If \code{NULL}, it defaults to 3.
#' @param clust_method Character. Clustering method applied to the first two
#'   FPCA scores to identify latent sub-structures. One of:
#'   \itemize{
#'     \item \code{"kmeans"} (default) — fast, effective for spherical
#'       clusters of similar size.
#'     \item \code{"hclust"} — agglomerative hierarchical clustering with
#'       Ward's D2 linkage; deterministic and robust to outliers.
#'     \item \code{"gmm"} — Gaussian Mixture Model via \code{mclust};
#'       handles elliptical clusters of varying shape and size, more realistic
#'       for biological data.
#'   }
#' @param factorNoise A numeric multiplier used to scale the adaptive noise injected into the simulated scores.
#' Higher values increase the simulated diversity (default is 1).
#' @param plot Logical. If \code{TRUE}, the function displays a 2 x 2 diagnostic dashboard
#' showing variance explanation, density zones, latent projections, and reconstructed spectra.
#'
#' @details
#' The function transforms the input mass spectrometry matrix into functional data objects using \code{MakeFPCAInputs}.
#' An FPCA is then performed to identify the main functional principal components with \code{AIC criterion}.
#'
#' Clustering is applied via the chosen \code{clust_method} on the first two
#' FPCA score dimensions to capture the underlying sub-structures.
#'
#' New scores are simulated cluster by cluster by adding an adaptive Gaussian noise proportional to the empirical variance of each component.
#'
#' Finally, spectra are reconstructed using the Karhunen-Loeve expansion formula:
#' \deqn{X_{s}(t) = \mu(t) + \sum_{k=1}^{K} \xi_{k, s} \phi_k(t)} where \eqn{\mu(t)} is the estimated mean function,
#' \eqn{\phi_k(t)} represent the eigenfunctions, and \eqn{\xi_{k, s}} are the newly simulated scores.
#' The reconstructed curves are then bounded to positive intensities via a flattening step and interpolated back onto the original
#' m/z grid.
#'
#' @return A numeric matrix with \code{n} rows and \code{p} columns containing the simulated spectra. Rows are named after their simulated cluster origin.
#'
#' @examples
#' \dontrun{
#' # Import 'spectra100' data simulate integrated into the package:
#' data(spectra100)
#'
#' # Default: K-means clustering
#' sim_kmeans <- SimulateSpectrum(data = spectra100, n = 50, k = 3)
#'
#' # Hierarchical clustering
#' sim_hclust <- SimulateSpectrum(data = spectra100, n = 50, k = 3, clust_method = "hclust")
#'
#' # Gaussian Mixture Model clustering
#' sim_gmm <- SimulateSpectrum(data = spectra100, n = 50, k = 3, clust_method = "gmm")
#'
#' # With diagnostic plots
#' SimulateSpectrum(data = spectra100, n = 60, k = 4, clust_method = "gmm", factorNoise = 1.2, plot = TRUE)
#' }
#'
#' @importFrom stats kmeans hclust cutree dist rnorm approx sd
#' @importFrom graphics par barplot lines text legend title points
#' @importFrom grDevices adjustcolor rainbow
#' @export

SimulateSpectrum <- function(data, n, k = NULL,
                             clust_method = c("kmeans", "hclust", "gmm"),
                             factorNoise = 1, plot = FALSE) {

  # -- Input validation
  clust_method <- match.arg(clust_method)
  if (is.null(k)) k <- 3L

  m_z            <- seq_len(ncol(data))
  colnames(data) <- m_z

  # -- 1. FPCA decomposition and clustering
  fpca_res      <- .fpcaClustering(data, k = k, n = n, clust_method = clust_method)
  res_fpca      <- fpca_res$res_fpca
  scores        <- fpca_res$scores
  cluster_ids   <- fpca_res$cluster_ids
  n_per_cluster <- fpca_res$n_per_cluster

  # -- 2. Score simulation with adaptive noise
  K_total         <- ncol(scores)
  sim_scores_list <- vector("list", k)

  for (i in seq_len(k)) {
    idx_base <- which(cluster_ids == i)
    pts_base <- scores[idx_base, , drop = FALSE]

    sds             <- apply(pts_base, 2, sd)
    sds[is.na(sds)] <- 0.0001

    n_sim <- n_per_cluster[i]
    if (n_sim > 0) {
      sim_idx <- sample(seq_len(nrow(pts_base)), n_sim, replace = TRUE)
      noise   <- matrix(rnorm(n_sim * K_total, mean = 0, sd = 1),
                        ncol = K_total)
      noise   <- sweep(noise, 2, sds * 0.1 * factorNoise, "*")
      sim_scores_list[[i]] <- pts_base[sim_idx, ] + noise
    }
  }

  all_sim_scores <- do.call(rbind, sim_scores_list)

  # -- 3. Spectral reconstruction
  simulated_matrix <- t(apply(
    all_sim_scores, 1, .reconstructSpectrum,
    mu          = res_fpca$mu,
    phi         = res_fpca$phi,
    targetGrid  = m_z,
    currentGrid = res_fpca$workGrid
  ))

  rownames(simulated_matrix) <- paste0("strain ", rep(seq_len(k), n_per_cluster))
  colnames(simulated_matrix) <- colnames(data)

  # -- 4. Diagnostic plots
  if (plot) {
    .simulateDiagPlots(
      res_fpca         = res_fpca,
      all_sim_scores   = all_sim_scores,
      simulated_matrix = simulated_matrix,
      m_z              = m_z,
      k                = k,
      n_per_cluster    = n_per_cluster
    )
  }

  return(simulated_matrix)
}
