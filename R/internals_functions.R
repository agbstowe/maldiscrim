# These functions are not exported and not visible in the package documentation.

# .fdaDecomposition

#' Functional decomposition of spectral data
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

#' Select the optimal number of PLS components via one-sigma heuristic
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

#' Project new spectral data into the functional space of a fitted model
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

#' Compute LDA posterior probabilities and alpha transparency values
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

#' Select a colour palette adapted to the number of groups
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


# .fpcaClustering

#' FPCA decomposition and clustering on spectral data
#'
#' Transforms a raw spectral matrix into functional data via
#' [fdapace::MakeFPCAInputs()] and [fdapace::FPCA()], then applies the chosen
#' clustering method on the first two functional principal component scores to
#' identify latent sub-structures.
#'
#' @param data A numeric matrix of spectral intensities (rows = samples,
#'   columns = m/z values). Column names must be set to the m/z indices before
#'   calling this function.
#' @param k Integer. Number of clusters. Must be >= 2.
#' @param n Integer. Total number of spectra to simulate. Used to compute
#'   `n_per_cluster` proportionally to observed cluster sizes.
#' @param clust_method Character. One of `"kmeans"` (default), `"hclust"`, or
#'   `"gmm"`. See Details.
#'
#' @details
#' The three supported clustering methods all operate on the first two FPCA
#' score dimensions:
#' \itemize{
#'   \item `"kmeans"` — K-means via [stats::kmeans()] with `nstart = 25` and
#'     a fixed seed for reproducibility. Fast and effective when clusters are
#'     roughly spherical and of similar size.
#'   \item `"hclust"` — Agglomerative hierarchical clustering via
#'     [stats::hclust()] with Ward's D2 linkage (`method = "ward.D2"`),
#'     followed by [stats::cutree()] at `k` groups. Deterministic (no seed
#'     needed) and robust to outliers.
#'   \item `"gmm"` — Gaussian Mixture Model via [mclust::Mclust()] with `G =
#'     k` components. Allows elliptical clusters of varying shape and size,
#'     which is more realistic for biological data. The covariance structure is
#'     selected automatically by BIC among standard `mclust` model families.
#' }
#'
#' @return A named list with the following elements:
#' \item{res_fpca}{The fitted FPCA object returned by [fdapace::FPCA()].}
#' \item{scores}{Numeric matrix of FPCA scores (rows = samples).}
#' \item{cluster_ids}{Integer vector of cluster assignments for each sample.}
#' \item{n_per_cluster}{Integer vector of target simulation counts per cluster,
#'   proportional to the observed cluster sizes. Adjusted so that the total
#'   equals `n`.}
#'
#' @keywords internal
.fpcaClustering <- function(data, k, n, clust_method = "kmeans") {

  m_z <- as.numeric(colnames(data))

  # -- FPCA
  input_data <- fdapace::MakeFPCAInputs(
    IDs  = rep(seq_len(nrow(data)), each = ncol(data)),
    tVec = rep(m_z, nrow(data)),
    yVec = as.vector(t(data))
  )

  res_fpca <- fdapace::FPCA(
    input_data$Ly, input_data$Lt,
    optns = list(dataType = "Dense", methodSelectK = "AIC")
  )
  scores      <- res_fpca$xiEst
  scores_2d   <- scores[, 1:2, drop = FALSE]

  # -- Clustering
  cluster_ids <- switch(clust_method,

                        "kmeans" = {
                          set.seed(42)
                          km <- kmeans(scores_2d, centers = k, nstart = 25)
                          km$cluster
                        },

                        "hclust" = {
                          d   <- dist(scores_2d)
                          hc  <- hclust(d, method = "ward.D2")
                          cutree(hc, k = k)
                        },

                        "gmm" = {
                          gmm <- mclust::Mclust(scores_2d, G = k, verbose = FALSE)
                          if (is.null(gmm)) {
                            stop(sprintf(
                              "GMM fitting failed for k = %d. Try a different k or clust_method.",
                              k
                            ))
                          }
                          gmm$classification
                        },

                        stop(sprintf(
                          "'clust_method' must be one of \"kmeans\", \"hclust\", or \"gmm\". Got: \"%s\".",
                          clust_method
                        ))
  )

  # -- Proportional allocation
  prop_clusters <- table(cluster_ids) / length(cluster_ids)
  n_per_cluster <- as.vector(round(prop_clusters * n))
  remainder     <- n - sum(n_per_cluster)
  if (remainder != 0) n_per_cluster[1] <- n_per_cluster[1] + remainder

  list(
    res_fpca      = res_fpca,
    scores        = scores,
    cluster_ids   = cluster_ids,
    n_per_cluster = n_per_cluster
  )
}

# .simulateDiagPlots ----------------------------------------------------------

#' Diagnostic plot dashboard for SimulateSpectrum (internal)
#'
#' Renders a 2x2 diagnostic dashboard summarising the FPCA simulation results:
#' \enumerate{
#'   \item Variance scree plot with cumulative variance overlay.
#'   \item KDE density zones in the FPCA latent space (real data).
#'   \item Simulated scores projected onto the KDE background.
#'   \item Reconstructed simulated spectra coloured by cluster.
#' }
#'
#' @param res_fpca The fitted FPCA object returned by [fdapace::FPCA()].
#' @param all_sim_scores Numeric matrix of simulated FPCA scores
#'   (rows = simulated samples).
#' @param simulated_matrix Numeric matrix of reconstructed simulated spectra
#'   (rows = simulated samples, columns = m/z values).
#' @param m_z Numeric vector of m/z indices.
#' @param k Integer. Number of clusters.
#' @param n_per_cluster Integer vector. Number of simulated samples per cluster.
#'
#' @return Called for its side effect (plots). Returns `NULL` invisibly.
#'
#' @keywords internal
.simulateDiagPlots <- function(res_fpca, all_sim_scores, simulated_matrix,
                               m_z, k, n_per_cluster) {

  K_total     <- ncol(res_fpca$xiEst)
  palette_sim <- rainbow(k)
  color_idx   <- rep(seq_len(k), n_per_cluster)
  strain_labels <- paste0("strain ", seq_len(k))

  # Save and restore par on exit — guards against fdapace resetting par
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  par(mfrow = c(2, 2), mar = c(4.5, 4.5, 2.5, 1.5))

  # A. Variance scree plot
  eigenvalues <- res_fpca$lambda
  var_exp     <- (eigenvalues / sum(eigenvalues)) * 100
  var_cum     <- cumsum(var_exp)

  bp <- barplot(
    var_exp,
    names.arg = paste0("FPC", seq_len(K_total)),
    col       = "steelblue",
    border    = "white",
    ylim      = c(0, 100),
    main      = "Variance Explained per Component",
    ylab      = "% Variance",
    cex.names = 0.7
  )
  lines(x = bp, y = var_cum, type = "b", pch = 19, col = "red", lwd = 1.5)
  text(x = bp, y = var_cum,
       labels = paste0(round(var_cum, 0), "%"),
       pos = 3, cex = 0.7, col = "red")

  # B. Density zones (KDE) — real data only
  # Snapshot par before fdapace call to restore mfrow/mar afterwards
  par_before_b <- par(no.readonly = TRUE)
  fdapace::CreateOutliersPlot(
    res_fpca,
    optns = list(fIndices = c(1, 2), variant = "KDE")
  )
  par(par_before_b)
  title(main = "Density Zones (Latent Space)", cex.main = 0.9)

  # C. Simulated scores projected onto KDE background --------------------------
  # Step 1: draw KDE background (real data)
  par_before_c <- par(no.readonly = TRUE)
  fdapace::CreateOutliersPlot(
    res_fpca,
    optns = list(fIndices = c(1, 2), variant = "KDE")
  )
  par(par_before_c)

  # Step 2: attenuate real data points by overlaying a semi-transparent white
  # rectangle — this dims the real points without removing the KDE contours
  usr <- par("usr")
  rect(usr[1], usr[3], usr[2], usr[4],
       col    = adjustcolor("white", alpha.f = 0.45),
       border = NA)

  # Step 3: overlay simulated points with high visibility
  points(
    all_sim_scores[, 1], all_sim_scores[, 2],
    pch = 3,
    col = adjustcolor(palette_sim[color_idx], alpha.f = 0.9),
    cex = 1.0
  )

  title(main = "Projection of Simulated Spectra", cex.main = 0.9)

  legend(
    "topright",
    legend = strain_labels,
    col    = palette_sim,
    pch    = 3,
    pt.cex = 1.0,
    cex    = 0.75,
    bty    = "n",                          # no legend box border
    title  = "Cluster"
  )

  # D. Reconstructed simulated spectra coloured by cluster ---------------------
  plot(m_z, simulated_matrix[1, ], type = "n",
       ylim = range(simulated_matrix),
       main = "Reconstructed Simulated Spectra",
       xlab = "m/z",
       ylab = "Intensity")

  for (i in seq_len(nrow(simulated_matrix))) {
    lines(m_z, simulated_matrix[i, ],
          col = adjustcolor(palette_sim[color_idx[i]], alpha.f = 0.3),
          lwd = 0.5)
  }

  legend(
    "topright",
    legend = strain_labels,
    col    = palette_sim,
    lty    = 1,
    lwd    = 1.5,
    cex    = 0.75,
    bty    = "n",
    title  = "Cluster"
  )

  invisible(NULL)
}
