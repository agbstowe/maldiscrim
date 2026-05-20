#' Visualisation functions for maldiscrim outputs
#'
#' @description
#' This group of functions provides a comprehensive suite of tools to visually explore, diagnose, and interpret the results of a fitted \code{fPLS_DA} model.
#' It includes classical 2D/3D score projections, UMAP manifold embeddings, variance scree plots, and spectral loading weight profiles.
#'
#' @param object A fitted model object of class inheriting from \code{"fPLS_DA"}.
#' @param comp A numeric vector specifying the projection axes or component to display.
#' @param show_proba Logical. If \code{TRUE}, sample points are colored by group, and their transparency is dynamically scaled based on their
#' maximum LDA posterior probability tier: 1.0 (\eqn{\ge 0.8}), 0.6 (\eqn{\ge 0.5}), or 0.2 (\eqn{< 0.5}).
#' @param ... Further arguments passed to or from other internal methods.
#'
#' @details
#' The score-based plotting functions (\code{plot_spectra_2d}, \code{plot_spectra_3d},
#' and \code{plot_spectra_umap}) assist in assessing cluster separation in different latent spaces.
#' When \code{show_proba = TRUE}, they act as a diagnostic layer to visually isolate samples near classification boundaries.
#'
#' Diagnostic plots (\code{plot_pls_var} and \code{plot_pls_weights}) help identify
#' the optimal number of components and track down specific mass-to-charge (\eqn{m/z}) channels.
#'
#' @examples
#' \dontrun{
#' # Assuming 'fit_model' is generated
#' # 1. 2D Score Projection
#' plot_spectra_2d(fit_model, comp = c(1, 2), show_proba = TRUE)
#'
#' # 2. Interactive 3D Plot
#' plot_spectra_3d(fit_model, comp = c(1, 2, 3))
#'
#' # 3. Variance Scree Plot
#' plot_pls_var(fit_model)
#'
#' # 4. Loading Profiles
#' plot_pls_weights(fit_model, comp = 1)
#'
#' # 5. Non-linear UMAP Mapping
#' plot_spectra_umap(fit_model, n_neighbors = 10, min_dist = 0.05)
#' }
#'
#' @importFrom stats predict
#' @name plot_functions
NULL

#' @rdname plot_functions
#' @section 3D Score Projection:
#' Generates an interactive tridimensional scatter plot of the PLS latent scores using \code{plotly}.
#' Axis titles automatically integrate the percentage of explained variance for each selected component.
#'
#' @export
plot_spectra_3d <- function(object, comp = c(1, 2, 3), show_proba = FALSE , ...) {
  scores_df <- as.data.frame(object$pls_model$scores[, comp])
  colnames(scores_df) <- c("X", "Y", "Z")
  scores_df$Souche <- as.factor(object$labels)

  pct <- round(object$pls_model$Xvar / object$pls_model$Xtotvar * 100, 2)

  # Gestion de la probabilité d'appartenance aux groupes.
  scores_df$Alpha <- 0.8 # Alpha par défaut
  if (show_proba) {
    raw_scores <- object$pls_model$scores[, 1:object$ncompOpt, drop = FALSE]
    clean_scores <- data.frame(raw_scores[, !object$const_idx, drop = FALSE])
    colnames(clean_scores) <- colnames(object$lda_model$means)

    lda_res <- stats::predict(object$lda_model, newdata = clean_scores)
    scores_df$Proba <- apply(lda_res$posterior, 1, max)
    scores_df$Alpha <- ifelse(scores_df$Proba >= 0.8, 1,
                       ifelse(scores_df$Proba >= 0.5, 0.6, 0.2))
  }

  # Calcul automatique des marges. On ajoute 10% d'espace sur chaque axe
  get_lims <- function(vec) {
    r <- range(vec)
    d <- diff(r)
    c(r[1] - 0.1 * d, r[2] + 0.1 * d)
  }

  xlims <- get_lims(scores_df$X)
  ylims <- get_lims(scores_df$Y)
  zlims <- get_lims(scores_df$Z)

  p <- plotly::plot_ly(scores_df, x = ~X, y = ~Y, z = ~Z, color = ~Souche,
                  type = 'scatter3d', mode = 'markers', colors = "Set1",
                  marker = list(size = 3, opacity = ~Alpha,...))

  p <- plotly::layout(p,scene = list(xaxis = list(title = paste0("Comp ", comp[1]," (", pct[comp[1]], "%)") , range = xlims),
                                yaxis = list(title = paste0("Comp ", comp[2]," (", pct[comp[2]], "%)") , range = ylims),
                                zaxis = list(title = paste0("Comp ", comp[3]," (", pct[comp[3]], "%)") , range = zlims)
                                ))

  return(p)
}


#' @rdname plot_functions
#' @section 2D Score Projection:
#' Generates a static two-dimensional scatter plot using \code{ggplot2}.
#'
#' @export
plot_spectra_2d <- function(object, comp = c(1, 2), show_proba = FALSE, ...) {
  # Construction du dataframe de base
  df <- data.frame(object$pls_model$scores[, comp])
  colnames(df) <- c("x", "y")
  df$Souche <- as.factor(object$labels)

  pct <- round(object$pls_model$Xvar / object$pls_model$Xtotvar * 100, 2)

  # On prend les min/max des composantes et on ajoute 10% de marge de chaque côté
  range_x <- range(df$x)
  range_y <- range(df$y)

  # Calcul des marges (expand)
  xlim <- c(range_x[1] - 0.1 * diff(range_x), range_x[2] + 0.1 * diff(range_x))
  ylim <- c(range_y[1] - 0.1 * diff(range_y ), range_y[2] + 0.1 * diff(range_y))

  # Initialisation du plot
  p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, color = Souche)) +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = paste("Comp", comp[1]," (", pct[comp[1]], "%)"), y = paste("Comp", comp[2], " (", pct[comp[2]], "%)")) +
    ggplot2::scale_color_brewer(palette = "Set1") +
    ggplot2::coord_cartesian(xlim = xlim, ylim = ylim)

  # Gestion de la probabilite si TRUE
  if (show_proba) {
    raw_scores <- object$pls_model$scores[, 1:object$ncompOpt, drop = FALSE]
    clean_scores <- data.frame(raw_scores[, !object$const_idx, drop = FALSE])
    colnames(clean_scores) <- colnames(object$lda_model$means)
    lda_res <- predict(object$lda_model, newdata = clean_scores)
    df$Proba <- apply(lda_res$posterior, 1, max)

    # On ajoute la couche avec alpha variable
    # p <- p + ggplot2::geom_point(data = df, ggplot2::aes(x = x, y = y, color = Souche, alpha = Proba), size = 2) +
    #  ggplot2::scale_alpha_continuous(limits = c(0, 1), range = c(0.2, 1))

    ########################
    # Création des paliers d'opacité demandés
    # On utilise cut() ou des ifelse pour transformer les probas en valeurs alpha discrètes
    df$Alpha <- ifelse(df$Proba >= 0.8, 1,
                            ifelse(df$Proba >= 0.5, 0.6, 0.2))

    # Ajout de la couche avec identité pour l'alpha
    p <- p + ggplot2::geom_point(data = df,
                                 ggplot2::aes(x = x, y = y, alpha = Alpha),
                                 size = 2) +
      ggplot2::scale_alpha_identity()

  } else {
    # Comportement par defaut (tout opaque)
    p <- p + ggplot2::geom_point(size = 2, alpha = 0.8)
  }

  return(p)
}



#' @rdname plot_functions
#' @section Explained Variance Scree Plot:
#' Computes and displays a bar chart of the individual variance percentage explained by each functional PLS component.
#'
#' @export
plot_pls_var <- function(object) {
  var_exp <- object$pls_model$Xvar / object$pls_model$Xtotvar * 100
  df <- data.frame(Comp = 1:length(var_exp), Var = var_exp)

  ggplot2::ggplot(df[1:object$ncompOpt,], ggplot2::aes(x = Comp, y = Var)) +
    ggplot2::geom_bar(stat = "identity", fill = "steelblue") +
    ggplot2::labs(title = "Variance expliquée par composante", y = "% Variance")
}


#' @rdname plot_functions
#' @section Spectral Loadings Profile:
#' Plots the loading weight profiles across variables (mass channels).
#' Useful for identifying specific \eqn{m/z} markers that highly influence class segregation.
#'
#' @export
plot_pls_weights <- function(object, comp = 1) {
  # Pour voir quelles masses (m/z) tirent la decision
  loading_values <- object$pls_model$loadings[, comp]
  df <- data.frame(m_z = 1:length(loading_values), Weight = loading_values)

  ggplot2::ggplot(df, ggplot2::aes(x = m_z, y = Weight)) +
    ggplot2::geom_line() +
    ggplot2::labs(title = paste("Influence des masses sur la Comp", comp))
}


#' @rdname plot_functions
#' @section UMAP Non-linear Embedding:
#' Projects the full matrix of high-dimensional PLS scores into a non-linear 2D layout via Uniform Manifold Approximation and Projection (UMAP).
#'
#' @param n_neighbors Integer. Number of nearest neighbors configuration for the UMAP algorithm.
#' @param min_dist Numeric. Controls how tightly UMAP packs points together in the layout.
#'
#' @export
plot_spectra_umap <- function(object, n_neighbors = 15, min_dist = 0.1, show_proba = FALSE, ...) {

  cf <- umap::umap.defaults
  cf$n_neighbors <- n_neighbors
  cf$min_dist <- min_dist

  umap_res <- umap::umap(as.matrix(unclass(object$pls_model$scores)),config = cf, ...)
  df <- data.frame(umap_res$layout)
  colnames(df) <- c("U1", "U2")
  df$Souche <- as.factor(object$labels)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = U1, y = U2, color = Souche))

  # Gestion de l'alpha via show_proba
  # alpha_val <- 1
  if (show_proba) {
    raw_scores <- object$pls_model$scores[, 1:object$ncompOpt, drop = FALSE]
    clean_scores <- data.frame(raw_scores[, !object$const_idx, drop = FALSE])
    colnames(clean_scores) <- colnames(object$lda_model$means)

    lda_res <- stats::predict(object$lda_model, newdata = clean_scores)
    df$Proba <- apply(lda_res$posterior, 1, max)
    df$Alpha <- ifelse(df$Proba >= 0.8, 1,
                       ifelse(df$Proba >= 0.5, 0.6, 0.2))

    # Création du plot avec transparence dynamique
    p <- p + ggplot2::geom_point(data = df, ggplot2::aes(alpha = Alpha), size = 2) +
      ggplot2::scale_alpha_identity()

  } else {
    # Plot standard
    p <- p + ggplot2::geom_point(size = 2, alpha = 0.8)
  }

  # Calcul automatique des marges (coord_cartesian)
  range_u1 <- range(df$U1)
  range_u2 <- range(df$U2)

  xlims <- c(range_u1[1] - 0.1 * diff(range_u1), range_u1[2] + 0.1 * diff(range_u1))
  ylims <- c(range_u2[1] - 0.1 * diff(range_u2), range_u2[2] + 0.1 * diff(range_u2))

  # Habillage du graphique
  p <- p +
    ggplot2::theme_light() +
    ggplot2::scale_color_brewer(palette = "Set1") +
    ggplot2::labs(title = "Projection UMAP des profils spectraux",
                  subtitle = paste0("Parameters : n_neighbors = ", n_neighbors, ", min_dist = ", min_dist),
                  x = "UMAP 1", y = "UMAP 2") +
    ggplot2::coord_cartesian(xlim = xlims, ylim = ylims)

  return(p)
}
