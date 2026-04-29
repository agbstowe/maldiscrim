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
#' @param comp Axes de projection des scores (par défaut 1, 2 et 3)
#' @export
plot_spectra_3d <- function(object, comp = c(1, 2, 3)) {
  scores_df <- as.data.frame(object$pls_model$scores[, comp])
  colnames(scores_df) <- c("X", "Y", "Z")
  scores_df$Souche <- as.factor(object$labels)

  p <- plotly::plot_ly(scores_df, x = ~X, y = ~Y, z = ~Z, color = ~Souche,
                  type = 'scatter3d', mode = 'markers')

  p <- plotly::layout(p,scene = list(xaxis = list(title = paste0("Comp ", comp[1])),
                                yaxis = list(title = paste0("Comp ", comp[2])),
                                zaxis = list(title = paste0("Comp ", comp[3]))
                                ))

  return(p)
}

#  PROJECTION 2D dans le plan des composantes de la PLS
#' @rdname plot_functions
#' @param object a fitted object of class inheriting from "fPLS_DA".
#' @param comp Axes de projection des scores (par défaut 1 et 2)
#' @param show_proba Booléen. Si TRUE, colorie par groupe et ajuste la transparence (alpha) selon la probabilité de prédiction.
#' @export
plot_spectra_2d <- function(object, comp = c(1, 2), show_proba = FALSE) {
  # Construction du dataframe de base
  df <- data.frame(object$pls_model$scores[, comp])
  colnames(df) <- c("x", "y")
  df$Souche <- as.factor(object$labels)

  # Initialisation du plot
  p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, color = Souche)) +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = paste("Comp", comp[1]), y = paste("Comp", comp[2])) +
    ggplot2::scale_color_brewer(palette = "Set1")

  # Gestion de la probabilité si TRUE
  if (show_proba) {
    # On prend uniquement le nombre de composantes utilisées pour l'entrainement (1:ncomp)
    raw_scores <- object$pls_model$scores[, 1:object$ncomp, drop = FALSE]
    # On applique le filtre des constantes (ceux qui ont été ignorés lors de l'entraînement du LDA)
    clean_scores <- data.frame(raw_scores[, !object$const_idx, drop = FALSE])
    # On s'assure que les noms de colonnes correspondent à ce que le LDA attend
    colnames(clean_scores) <- colnames(object$lda_model$means)
    # On demande au LDA de nous donner les probas sur ces scores
    lda_res <- predict(object$lda_model, newdata = clean_scores)
    # On calcule la probabilité max pour chaque spectre
    df$Proba <- apply(lda_res$posterior, 1, max)

    # On ajoute la couche avec alpha variable
    p <- p + ggplot2::geom_point(data = df, ggplot2::aes(x = x, y = y, color = Souche, alpha = Proba), size = 2) +
      ggplot2::scale_alpha_continuous(limits = c(0, 1), range = c(0.2, 1))
  } else {
    # Comportement par défaut (tout opaque)
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
  # Pour voir quelles masses (m/z) tirent la décision
  loading_values <- object$pls_model$loadings[, comp]
  df <- data.frame(m_z = 1:length(loading_values), Weight = loading_values)

  ggplot2::ggplot(df, ggplot2::aes(x = m_z, y = Weight)) +
    ggplot2::geom_line() +
    ggplot2::labs(title = paste("Influence des masses sur la Comp", comp))
}

#  UMAP sur les scores PLS du modele de discrimination
#' @rdname plot_functions
#' @export
plot_spectra_umap <- function(object, n_neighbors = 15, min_dist = 0.1) {
  # On utilise les scores PLS comme base pour l'UMAP (plus robuste que les spectres bruts)
  umap_res <- umap::umap(as.matrix(unclass(object$pls_model$scores)), n_neighbors = n_neighbors, min_dist = min_dist)
  df <- data.frame(umap_res$layout)
  colnames(df) <- c("U1", "U2")
  df$Souche <- as.factor(object$labels)

  ggplot2::ggplot(df, ggplot2::aes(U1, U2, color = Souche)) +
    ggplot2::geom_point(size = 1, alpha = 0.8) +
    ggplot2::theme_light() +
    ggplot2::labs(title = "Projection UMAP des profils spectraux")
}
