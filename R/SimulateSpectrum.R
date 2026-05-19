#' Simulate spectrum using Functionnal Principal Components Analysis (FPCA)
#'
#' @description
#' Simulates a set of $n$ standard spectra based on real or simulated data using Functional Principal Components Analysis (FPCA) and adaptive noise injection.
#'
#' @param data A numeric matrix of Mass Spectrometry (MSP) data, where rows represent samples and columns represent m/z variables.
#' @param n Integer. Total number of spectra to simulate in the output matrix.
#' @param k Integer. Number of clusters/strains to identify in the latent space. If \code{NULL}, it defaults to 3.
#' @param factor_noise A numeric multiplier used to scale the adaptive noise injected into the simulated scores.
#' Higher values increase the simulated diversity (default is 1).
#' @param plot Logical. If \code{TRUE}, the function displays a $2 x 2$ diagnostic dashboard
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
#' @return A numeric matrix with $n$ rows and $p$ columns containing the simulated spectra. Rows are named after their simulated cluster origin.
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
#'                                             factor_noise = 1.2, plot = TRUE)
#' dim(simulated_data_plots)
#' }
#'
#' @importFrom stats kmeans rnorm approx sd
#' @importFrom graphics par barplot lines text legend title points
#' @importFrom grDevices adjustcolor rainbow
#' @export

SimulateSpectrum <- function(data, n, k = NULL, factor_noise = 1, plot = FALSE) {

  # 1. Préparation et FPCA
  m_z = 1:ncol(data)
  colnames(data) = m_z
  input_data <- fdapace::MakeFPCAInputs(IDs = rep(1:nrow(data), each = ncol(data)),
                                        tVec = rep(as.numeric(m_z), nrow(data)),
                                        yVec = as.vector(t(data)))

  res_fpca <- fdapace::FPCA(input_data$Ly, input_data$Lt, optns = list(dataType = 'Dense', methodSelectK = "AIC"))
  scores <- res_fpca$xiEst

  # 2. Clustering
  if (is.null(k)) { k <- 3 }
  set.seed(42)
  km <- kmeans(scores[, 1:2], centers = k, nstart = 25)
  cluster_ids <- km$cluster

  # 3. Calcul des effectifs cibles (Inchangé)
  prop_clusters <- table(cluster_ids) / length(cluster_ids)
  n_per_cluster <- as.vector(round(prop_clusters * n))
  diff <- n - sum(n_per_cluster)
  if (diff != 0) n_per_cluster[1] <- n_per_cluster[1] + diff

  # 4. Simulation avec BRUIT ADAPTATIF sur K composantes. On intègre toutes les composantes pour la simulation
  K_total <- ncol(scores)
  sim_scores_list <- list()

  for (i in 1:k) {
    idx_base <- which(cluster_ids == i)
    pts_base <- scores[idx_base, , drop = FALSE]

    # Calcul de l'écart-type pour CHAQUE composante (1 à K)
    # On simule la variabilité réelle de chaque axe pour ce cluster
    sds <- apply(pts_base, 2, sd)
    sds[is.na(sds)] <- 0.0001 # Sécurité

    n_sim <- n_per_cluster[i]
    if (n_sim > 0) {
      sim_idx <- sample(1:nrow(pts_base), n_sim, replace = TRUE)
      # Génération d'une matrice de bruit (n_sim x K_total)
      noise <- matrix(rnorm(n_sim * K_total, mean = 0, sd = 1), ncol = K_total)
      # On multiplie chaque colonne de bruit par le SD correspondant (adaptatif)
      # sweep est plus propre que cbind ici pour gérer K colonnes
      noise <- sweep(noise, 2, sds * 0.1 * factor_noise, "*")

      sim_scores_list[[i]] <- pts_base[sim_idx, ] + noise
    }
  }
  all_sim_scores <- do.call(rbind, sim_scores_list)

  # 5. Reconstruction avec redressement (K-composantes)
  mu <- res_fpca$mu
  phi <- res_fpca$phi
  workGrid <- res_fpca$workGrid

  reconstruct_spectrum <- function(s, mu, phi, targetGrid, currentGrid) {
    # phi %*% s calcule automatiquement la somme Mu + Score1*Phi1 + ... + ScoreK*PhiK
    y_grid <- as.vector(mu + (phi %*% s))

    # Redressement physique
    y_grid <- pmax(y_grid, 0)

    # Interpolation sur les Daltons d'origine (ex: 20 664 points)
    y_final <- approx(x = currentGrid, y = y_grid, xout = targetGrid)$y
    return(y_final)
  }

  simulated_matrix <- t(apply(all_sim_scores, 1, reconstruct_spectrum,
                              mu = mu, phi = phi,
                              targetGrid = m_z,
                              currentGrid = workGrid))

  # Finalisation des noms
  # On utilise rep(1:k, n_per_cluster) pour que les labels collent aux blocs de la liste
  rownames(simulated_matrix) <- paste("souche ", rep(1:k, n_per_cluster))
  colnames(simulated_matrix) <- colnames(data)


  # 6. GRAPHIQUES DE DIAGNOSTIC
  if (plot) {
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par))

    # Configuration : 2x2 graphiques avec des marges standards
    par(mfrow = c(2, 2), mar = c(4.5, 4.5, 2.5, 1.5))

    # A. Eboulis des variances (Inchangé mais sans légende)
    eigenvalues <- res_fpca$lambda
    var_exp <- (eigenvalues / sum(eigenvalues)) * 100
    var_cum <- cumsum(var_exp)
    bp <- barplot(var_exp, names.arg = paste0("FPC", 1:K_total),
                  col = "steelblue", border = "white", ylim = c(0, 100),
                  main = "Variance Expliquee par Composante", ylab = "% Variance", cex.names = 0.7)
    lines(x = bp, y = var_cum, type = "b", pch = 19, col = "red", lwd = 1.5)
    text(x = bp, y = var_cum, labels = paste0(round(var_cum, 0), "%"), pos = 3, cex = 0.7, col = "red")
    # Legend supprimee

    # B. Zones de Densite (KDE) - Seulement le réel pour le fond
    fdapace::CreateOutliersPlot(res_fpca, optns = list(fIndices = c(1,2), variant = "KDE"))
    title(main = "Zones de Densite (Espace Latent)", cex.main = 0.9)
    # Legend de fdapace supprimee

    # C. Projection des spectres simules SUR le fond createoutliersplot
    # 1. On trace le fond KDE
    fdapace::CreateOutliersPlot(res_fpca, optns = list(fIndices = c(1,2), variant = "KDE"))

    # 2. On superpose les points simulés en forçant le tracé sur le graphique existant
    # Nous masquons les points réels en utilisant le même fond KDE

    # On définit les couleurs (Rainbow est plus automatique pour le package)
    palette_sim <- rainbow(k)

    # 3. Superposition des points simulés (+) sans points réels
    points(all_sim_scores[,1], all_sim_scores[,2],
           pch = 3, col = palette_sim[rep(1:k, n_per_cluster)],
           cex = 0.8)

    title(main = "Projection des spectres simules", cex.main = 0.9)
    # Legend supprimee

    # D. Spectres Reconstitues (Colories par souche, sans légende)
    plot(m_z, simulated_matrix[1,], type = "n",
         ylim = range(simulated_matrix),
         main = "Spectres Simules Reconstitues", xlab = "m/z", ylab = "Intensite")

    for(i in 1:nrow(simulated_matrix)) {
      #On utilise l'index numerique pour la couleur
      color_idx <- rep(1:k, n_per_cluster)
      lines(m_z, simulated_matrix[i,],
            col = adjustcolor(palette_sim[color_idx[i]], alpha.f = 0.3), lwd = 0.5)
    }

  }

  return(simulated_matrix)

}
