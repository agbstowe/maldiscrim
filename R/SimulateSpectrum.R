#' Simulate spectrum using Functionnal Principal Components Analysis (FPCA)
#'
#' @description
#' Simulates a set of n standard spectra based on real or simulated data using Functional Principal Components Analysis (FPCA) and adaptive noise injection.
#'
#' @param data A numeric matrix of Mass Spectrometry (MSP) data, where rows represent samples and columns represent m/z variables.
#' @param n Integer. Total number of spectra to simulate in the output matrix.
#' @param k Integer. Number of clusters/strains to identify in the latent space. If \code{NULL}, it defaults to 3.
#' @param factorNoise A numeric multiplier used to scale the adaptive noise injected into the simulated scores.
#' Higher values increase the simulated diversity (default is 1).
#' @param plot Logical. If \code{TRUE}, the function displays a 2 x 2 diagnostic dashboard
#' showing variance explanation, density zones, latent projections, and reconstructed spectra.
#'
#' @details
#' The function transforms the input mass spectrometry matrix into functional data objects using \code{MakeFPCAInputs}.
#' An FPCA is then performed to identify the main functional principal components with \code{AIC criterion}.
#' Clustering is applied via K-means on the first two principal component scores to capture the underlying sub-structures.
#'
#' New scores are simulated cluster by cluster by adding an adaptive Gaussian noise proportional to the empirical variance of each component.
#'
#' Finally, spectra are reconstructed using the Karhunen-Loeve expansion formula:
#' \deqn{X_{s}(t) = \mu(t) + \sum_{k=1}^{K} \xi_{k, s} \phi_k(t)} where \eqn{\mu(t)} is the estimated mean function,
#' \eqn{\phi_k(t)} represent the eigenfunctions, and \eqn{\xi_{k, s}} are the newly simulated scores.
#' The reconstructed curves are then bounded to positive intensities via a flattening step and interpolated back onto the original
#' m/z grid.
#'
#' @return A numeric matrix with n rows and p columns containing the simulated spectra. Rows are named after their simulated cluster origin.
#'
#' @examples
#' \dontrun{
#' # Import 'spectra100' data simulate integrated into the package:
#' data(spectra100)
#'
#' # Run a simulation of 50 spectra without plotting
#' simulated_data <- SimulateSpectrum(data = spectra100, n = 50, k = 3, plot = FALSE)
#' dim(simulated_data)
#'
#' # Run a simulation with diagnostic plots enabled
#' simulated_data_plots <- SimulateSpectrum(data = spectra100, n = 60, k = 4,
#'                                             factorNoise = 1.2, plot = TRUE)
#' dim(simulated_data_plots)
#' }
#'
#' @importFrom stats kmeans rnorm approx sd
#' @importFrom graphics par barplot lines text legend title points
#' @importFrom grDevices adjustcolor rainbow
#' @export

SimulateSpectrum <- function(data, n, k = NULL, factorNoise = 1, plot = FALSE) {

  # FPCA
  m_z <- seq_len(ncol(data))
  colnames(data) <- m_z

  input_data <- fdapace::MakeFPCAInputs(IDs = rep(seq_len(nrow(data)), each = ncol(data)),
                                        tVec = rep(as.numeric(m_z), nrow(data)),
                                        yVec = as.vector(t(data)))

  res_fpca <- fdapace::FPCA(input_data$Ly, input_data$Lt, optns = list(dataType = 'Dense', methodSelectK = "AIC"))
  scores <- res_fpca$xiEst

  # Clustering
  if (is.null(k)) { k <- 3L }
  set.seed(42)
  km          <- kmeans(scores[, 1:2], centers = k, nstart = 25)
  cluster_ids <- km$cluster

  prop_clusters <- table(cluster_ids) / length(cluster_ids)
  n_per_cluster <- as.vector(round(prop_clusters * n))
  remainder     <- n - sum(n_per_cluster)
  if (remainder != 0) n_per_cluster[1] <- n_per_cluster[1] + remainder

  # Simulation with adaptative noise
  K_total <- ncol(scores)
  sim_scores_list <- vector("list", k)

  for (i in seq_len(k)) {
    idx_base <- which(cluster_ids == i)
    pts_base <- scores[idx_base, , drop = FALSE]

    sds             <- apply(pts_base, 2, sd)
    sds[is.na(sds)] <- 0.0001

    n_sim <- n_per_cluster[i]
    if (n_sim > 0) {
      sim_idx <- sample(1:nrow(pts_base), n_sim, replace = TRUE)
      # noise matrix (n_sim x K_total)
      noise <- matrix(rnorm(n_sim * K_total, mean = 0, sd = 1), ncol = K_total)
      noise <- sweep(noise, 2, sds * 0.1 * factorNoise, "*")

      sim_scores_list[[i]] <- pts_base[sim_idx, ] + noise
    }
  }
  all_sim_scores <- do.call(rbind, sim_scores_list)

  # Spectral reconstruction with physical bounding
  mu       <- res_fpca$mu
  phi      <- res_fpca$phi
  workGrid <- res_fpca$workGrid

  reconstruct_spectrum <- function(s, mu, phi, targetGrid, currentGrid) {
    # phi %*% s computes mu + score_1 * phi_1 + ... + score_K * phi_K
    y_grid <- as.vector(mu + (phi %*% s))
    # Physical bounding: force non-negative intensities
    y_grid <- pmax(y_grid, 0)

    # Interpolate back onto the original m/z grid
    y_final <- approx(x = currentGrid, y = y_grid, xout = targetGrid)$y
    return(y_final)
  }

  simulated_matrix <- t(apply(all_sim_scores, 1, reconstruct_spectrum,
                              mu = mu, phi = phi,
                              targetGrid = m_z,
                              currentGrid = workGrid))

  rownames(simulated_matrix) <- paste("strain ", rep(seq_len(k), n_per_cluster))
  colnames(simulated_matrix) <- colnames(data)


  # Diagnostic plots
  if (plot) {
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par))
    par(mfrow = c(2, 2), mar = c(4.5, 4.5, 2.5, 1.5))

    # A. Variance scree plot
    eigenvalues <- res_fpca$lambda
    var_exp     <- (eigenvalues / sum(eigenvalues)) * 100
    var_cum     <- cumsum(var_exp)
    bp          <- barplot(var_exp, names.arg = paste0("FPC", seq_len(K_total)),
                  col = "steelblue", border = "white", ylim = c(0, 100),
                  main      = "Variance Explained per Component",
                  ylab      = "% Variance", cex.names = 0.7)
    lines(x = bp, y = var_cum, type = "b", pch = 19, col = "red", lwd = 1.5)
    text(x = bp, y = var_cum, labels = paste0(round(var_cum, 0), "%"), pos = 3, cex = 0.7, col = "red")

    # B. Density zones (KDE) — real data
    fdapace::CreateOutliersPlot(res_fpca, optns = list(fIndices = c(1,2), variant = "KDE"))
    title(main = "Density Zones (Latent Space)", cex.main = 0.9)

    # C. Simulated scores projected onto the KDE background
    fdapace::CreateOutliersPlot(res_fpca, optns = list(fIndices = c(1,2), variant = "KDE"))

    palette_sim <- rainbow(k)
    points(all_sim_scores[,1], all_sim_scores[,2],
           pch = 3, col = palette_sim[rep(seq_len(k), n_per_cluster)],
           cex = 0.8)

    title(main = "Projection of Simulated Spectra", cex.main = 0.9)

    # D. Reconstructed simulated spectra (coloured by cluster)
    plot(m_z, simulated_matrix[1,], type = "n",
         ylim = range(simulated_matrix),
         main = "Reconstructed Simulated Spectra",
         xlab = "m/z", ylab = "Intensity")

    color_idx <- rep(seq_len(k), n_per_cluster)
    for (i in seq_len(nrow(simulated_matrix))) {
      lines(m_z, simulated_matrix[i, ],
            col = adjustcolor(palette_sim[color_idx[i]], alpha.f = 0.3),
            lwd = 0.5)
    }

  }

  return(simulated_matrix)

}
