#' Visualisations des resultats de l'analyse discriminante
#' PROJECTION 3D dans l'espace des composantes de la PLS
#' PROJECTION 2D dans le plan des composantes de la PLS
#' Visualisation graphique de la VARIANCE EXPLIQUEE des composantes PLS (Scree Plot)
#  UMAP sur les scores PLS du modele de discrimination
#' @param object Un objet de classe 'spectra_model'
#' @param ... further arguments passed to or from other methods.
#' @name plot_functions
NULL

# PROJECTION 3D dans l'espace des composantes de la PLS
#' @rdname plot_functions
#' @param comp Axes de projection des scores (par defaut 1, 2 et 3)
#' @param show_proba Booleen. Si TRUE, colorie par groupe et ajuste la transparence (alpha) selon la probabilite de prediction.
#' @param ... further arguments passed to or from other methods like ggplot graphic topics.
#' @export
plot_spectra_3d <- function(object, comp = c(1, 2, 3), show_proba = FALSE , ...) {
  scores_df <- as.data.frame(object$pls_model$scores[, comp])
  colnames(scores_df) <- c("X", "Y", "Z")
  scores_df$Souche <- as.factor(object$labels)

  # Gestion de la probabilité d'appartenance aux groupes.
  alpha_val <- 1 # Alpha par défaut
  if (show_proba) {
    raw_scores <- object$pls_model$scores[, 1:object$ncomp, drop = FALSE]
    clean_scores <- data.frame(raw_scores[, !object$const_idx, drop = FALSE])
    colnames(clean_scores) <- colnames(object$lda_model$means)

    lda_res <- stats::predict(object$lda_model, newdata = clean_scores)
    scores_df$Proba <- apply(lda_res$posterior, 1, max)
    alpha_val <- scores_df$Proba # L'alpha devient dynamique
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
                  marker = list(size = 2,...))

  p <- plotly::layout(p,scene = list(xaxis = list(title = paste0("Comp ", comp[1]), range = xlims),
                                yaxis = list(title = paste0("Comp ", comp[2]), range = ylims),
                                zaxis = list(title = paste0("Comp ", comp[3]), range = zlims)
                                ))

  return(p)
}



#  PROJECTION 2D dans le plan des composantes de la PLS
#' @rdname plot_functions
#' @param object a fitted object of class inheriting from "fPLS_DA".
#' @param comp Axes de projection des scores (par défaut 1 et 2)
#' @param show_proba Booleen. Si TRUE, colorie par groupe et ajuste la transparence (alpha) selon la probabilite de prediction.
#' @param ... further arguments passed to or from other methods like ggplot graphic topics.
#' @export
plot_spectra_2d <- function(object, comp = c(1, 2), show_proba = FALSE, ...) {
  # Construction du dataframe de base
  df <- data.frame(object$pls_model$scores[, comp])
  colnames(df) <- c("x", "y")
  df$Souche <- as.factor(object$labels)

  # On prend les min/max des composantes et on ajoute 10% de marge de chaque côté
  range_x <- range(df$x)
  range_y <- range(df$y)

  # Calcul des marges (expand)
  xlim <- c(range_x[1] - 0.1 * diff(range_x), range_x[2] + 0.1 * diff(range_x))
  ylim <- c(range_y[1] - 0.1 * diff(range_y), range_y[2] + 0.1 * diff(range_y))

  # Initialisation du plot
  p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, color = Souche)) +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = paste("Comp", comp[1]), y = paste("Comp", comp[2])) +
    ggplot2::scale_color_brewer(palette = "Set1") +
    # ON COORD_CARTESIAN avec les marges calculées pour l'affichage par défaut
    ggplot2::coord_cartesian(xlim = xlim, ylim = ylim)

  # Gestion de la probabilite si TRUE
  if (show_proba) {
    # On prend uniquement le nombre de composantes utilisees pour l'entrainement (1:ncomp)
    raw_scores <- object$pls_model$scores[, 1:object$ncomp, drop = FALSE]
    # On applique le filtre des constantes (ceux qui ont ete ignores lors de l'entrainement du LDA)
    clean_scores <- data.frame(raw_scores[, !object$const_idx, drop = FALSE])
    # On s'assure que les noms de colonnes correspondent a ce que le LDA attend
    colnames(clean_scores) <- colnames(object$lda_model$means)
    # On demande au LDA de nous donner les probas sur ces scores
    lda_res <- predict(object$lda_model, newdata = clean_scores)
    # On calcule la probabilite max pour chaque spectre
    df$Proba <- apply(lda_res$posterior, 1, max)

    # On ajoute la couche avec alpha variable
    p <- p + ggplot2::geom_point(data = df, ggplot2::aes(x = x, y = y, color = Souche, alpha = Proba), size = 2) +
      ggplot2::scale_alpha_continuous(limits = c(0, 1), range = c(0.2, 1))
  } else {
    # Comportement par defaut (tout opaque)
    p <- p + ggplot2::geom_point(size = 2, alpha = 1)
  }

  return(p)
}



#  Visualisation graphique de la VARIANCE EXPLIQUEE des composantes PLS (Scree Plot)
#' @rdname plot_functions
#' @export
plot_pls_var <- function(object) {
  var_exp <- object$pls_model$Xvar / object$pls_model$Xtotvar * 100
  df <- data.frame(Comp = 1:length(var_exp), Var = var_exp)

  ggplot2::ggplot(df[1:object$ncomp,], ggplot2::aes(x = Comp, y = Var)) +
    ggplot2::geom_bar(stat = "identity", fill = "steelblue") +
    ggplot2::labs(title = "Variance expliquée par composante", y = "% Variance")
}



#  Visualisation de la representativite des variables dans la construction des composantes de la PLS (Loadings/Weights)
#' @rdname plot_functions
#' @export
plot_pls_weights <- function(object, comp = 1) {
  # Pour voir quelles masses (m/z) tirent la decision
  loading_values <- object$pls_model$loadings[, comp]
  df <- data.frame(m_z = 1:length(loading_values), Weight = loading_values)

  ggplot2::ggplot(df, ggplot2::aes(x = m_z, y = Weight)) +
    ggplot2::geom_line() +
    ggplot2::labs(title = paste("Influence des masses sur la Comp", comp))
}



#  UMAP sur les scores PLS du modele de discrimination
#' @rdname plot_functions
#' @param object a fitted object of class inheriting from "fPLS_DA".
#' @inheritParams umap::umap
#' @param show_proba Booleen. Si TRUE, colorie par groupe et ajuste la transparence (alpha) selon la probabilite de prediction.
#' @param ... further arguments passed to or from other methods.
#' @export
plot_spectra_umap <- function(object, n_neighbors = 15, min_dist = 0.1, show_proba = FALSE, ...) {
  # On utilise les scores PLS comme base pour l'UMAP (plus robuste que les spectres bruts)
  umap_res <- umap::umap(as.matrix(unclass(object$pls_model$scores)), n_neighbors = n_neighbors, min_dist = min_dist, ...)
  df <- data.frame(umap_res$layout)
  colnames(df) <- c("U1", "U2")
  df$Souche <- as.factor(object$labels)

  # Gestion de l'alpha via show_proba
  # alpha_val <- 1
  if (show_proba) {
    raw_scores <- object$pls_model$scores[, 1:object$ncomp, drop = FALSE]
    clean_scores <- data.frame(raw_scores[, !object$const_idx, drop = FALSE])
    colnames(clean_scores) <- colnames(object$lda_model$means)

    lda_res <- stats::predict(object$lda_model, newdata = clean_scores)
    df$Proba <- apply(lda_res$posterior, 1, max)

    # Création du plot avec transparence dynamique
    p <- ggplot2::ggplot(df, ggplot2::aes(U1, U2, color = Souche)) +
      ggplot2::geom_point(ggplot2::aes(alpha = Proba), size = 2) +
      ggplot2::scale_alpha_continuous(limits = c(0, 1), range = c(0.2, 1))
  } else {
    # Plot standard
    p <- ggplot2::ggplot(df, ggplot2::aes(U1, U2, color = Souche)) +
      ggplot2::geom_point(size = 2, alpha = 1)
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
                  x = "U1", y = "U2") +
    ggplot2::coord_cartesian(xlim = xlims, ylim = ylims)

  return(p)
}

