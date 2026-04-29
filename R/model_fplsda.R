#' Fit Partiel Least Square Regression and Discriminant Analysis on functionnal data (fPLS_DA)
#'
#' @param data Matrice de donnees MSP transformees pour l'entrainement du modele.
#' @param method Méthode de décomposition en bases fonctionnelles : "bsplines" ou "wavelets".
#' @param nbasis an integer variable specifying the number of basis functions (1050 par défaut pour bsplines).
#' @param ncomp Nombre de composantes PLS à retenir (par défaut 20).
#' @param argvals description
#' @param rangeval description
#' @param n_folds Nombre de plis pour la validation croisée (par défaut 5).
#' @inheritParams wavelets::dwt
#' @export
fit.fPLS_DA <- function(data,
                          method = c("bsplines", "wavelets"),
                          nbasis = 1050,
                          ncomp = NULL,
                          argvals = 2001:22664,
                          rangeval = c(1999,22664),
                          n_folds = 5,
                          filter = "la14",
                          boundary = "periodic",
                          n.levels = NULL,
                          level = 3) {

  method <- match.arg(method)
  groups <- as.factor(rownames(data))
  # if (is.null(argvals)){
  #   argvals <- 2001:22664  ## Grille de projection de tous les spectres. Tous les spectres doivent être calibrés sur cet intervalle
  # } else {
  #   argvals = argvals
  # }


  # 1. Décomposition selon la méthode choisie
  if (method == "bsplines") {
    message("Utilisation de la décomposition fonctionnelle en B-splines.")
    basis_obj <- fda::create.bspline.basis(rangeval = rangeval, nbasis = nbasis)
    fd_obj <- fda::Data2fd(argvals = argvals, y = t(data), basisobj = basis_obj)
    coef_matrix <- t(fd_obj$coefs)
    }
  else if (method == "wavelets") {
    message("Utilisation de la décomposition fonctionnelle en ondelettes.")

    max_levels <- floor(log2(ncol(data)))

    if (is.null(n.levels)) {
      # On prend la puissance de 2 la plus proche de la longueur du signal
      # log2(20664) donne environ 14.3, donc floor() nous donne 14
      n.levels <- max_levels
    } else if (n.levels > max_levels) {
      warning(paste("n.levels trop grand. Réglé automatiquement sur", max_levels))
      n.levels <- max_levels
    }

    # 3. Validation de la pertinence du 'level' (Garde-fou scientifique)
    # On s'assure qu'au niveau choisi, il reste assez de coefficients (ici > 20 pour éviter que la PLS ne plante)
    # Nombre théorique de coeffs ~= Longueur / 2^level
    n_coef_theorique <- floor(ncol(data) / (2^level))

    if (n_coef_theorique < 20) {
      new_level <- max(1, floor(log2(ncol(data) / 20)))
      warning(paste("Le level choisi ", level, " est trop élevé pour satisfaire la PLS. Ajustement automatique à : ", new_level))
      level <- new_level
    }

    if (level > n.levels) {
      warning(paste("Le level choisi ", level, " supérieur au nombre total de niveaux. Ajustement à ", n.levels))
      level <- n.levels
    }

    # On applique la DWT sur chaque spectre (chaque ligne de data)
    # On utilise 'boundary="periodic"' pour éviter les effets de bord au début/fin du spectre
    # On utilise apply pour créer une liste contenant les objets dwt
    wt_list <- apply(data, 1, function(x) {
      wavelets::dwt(x, filter = filter, n.levels = n.levels, boundary = boundary)
    })

    # Extraction des coefficients pour la PLS
    coef_matrix <- t(sapply(wt_list, function(wt) wt@W[[level]]))
    wt_output <- wt_list # On stocke la liste complète
  }

  # 2. PLS sur ncomp composantes et Gestion de la variance nulle
  y_dummy <- fastDummies::dummy_cols(groups, remove_selected_columns = TRUE)

  max_search <- min(nrow(data) - 1, 20) # On cherche jusqu'à 20 composantes max pour la pls pour éviter un éventuelle sur apprentissage

  pls_mod <- pls::plsr(as.matrix(y_dummy) ~ coef_matrix, ncomp = max_search, method = "kernelpls",
                       validation = "CV", segments = n_folds)

  # Optimisation automatique du nombre de composante
  if (!is.null(ncomp)) {
    # L'utilisateur a fourni une valeur pour ncomp et on borne la valeur par max_search pour éviter les plantages.
    ncomp <- min(ncomp, max_search)
    message("Nombre de composantes optimisé automatiquement à ", ncomp, " pour la PLS")

  } else {
    # Calcul automatique du ncomp optimal via la méthode onesigma de la fonction selectNcomp
    rmsep_data <- pls::RMSEP(pls_mod, estimate = "CV")$val
    n_responses <- dim(rmsep_data)[1]

    ncomp_vec <- sapply(1:n_responses, function(i) {
      errors <- drop(rmsep_data[i, 1, -1])
      min_idx <- which.min(errors)
      min_rmse <- errors[min_idx]
      # Calcul du seuil onesigma
      target <- min_rmse + sd(errors) / sqrt(length(errors))
      best_n <- which(errors <= target)[1]
      return(best_n)
    })

    ncomp <- max(ncomp_vec)
    message("Nombre de composantes optimisé automatiquement à ", ncomp, " pour la PLS")
  }

  # Calcul du filtrage de variance sur les scores
  scores <- pls_mod$scores[, 1:ncomp]
  group_means <- aggregate(scores, list(groups), mean)[,-1]
  f1 <- apply(scores - as.matrix(group_means[groups, ]), 2, sd)
  const_idx <- f1 < 0.0001

  # 3. Modèle LDA
  scores_clean <- scores[, !const_idx]
  lda_mod <- MASS::lda(x = scores_clean, grouping = groups)

  # 4. Construction de l'objet S3
  model_obj <- list(
    pls_model = pls_mod,
    lda_model = lda_mod,
    method = method,
    basis = if(method == "bsplines") basis_obj else NULL,
    wavelets = if(method == "wavelets") wt_output else NULL,
    ncomp = ncomp,
    const_idx = const_idx,
    levels = levels(groups),
    labels = groups,
    argvals = if(method == "bsplines") argvals else NULL,
    filter = if(method == "wavelets") filter else NULL,
    n_levels = if(method == "wavelets") n.levels else NULL,
    level = if(method == "wavelets") level else NULL,
    nbasis = if(method == "bsplines") nbasis else NULL,
    boundary = if(method == "wavelets") boundary else NULL
  )

  class(model_obj) <- "fPLS_DA"
  return(model_obj)
}

#' Predict Method for fPLS_DA fits
#' @param object a fitted object of class inheriting from "fPLS_DA".
#' @param newdata an optionally data frame in which to look for variables with which to predict.
#' If omitted, the fitted values are used.   ###################" Revoir le newdata et le rendre optionnel
#' @param threshold Seuil de validation du type de la souche
#' @param ... further arguments passed to or from other methods.
#'
#' @export
predict.fPLS_DA <- function(object, newdata, threshold = NULL) {

  # VALIDATION DÉFENSIF
  if (!is.matrix(newdata) && !is.data.frame(newdata)) {
    stop("newdata doit être une matrice ou un data.frame.")
  }
  # # Vérification cruciale : la taille du spectre. On compare avec la longueur de argvals
  # if (ncol(newdata) != length(object$argvals)) {
  #   stop(paste("Erreur : newdata possède", ncol(newdata),
  #              "colonnes, mais le modèle a été entraîné avec", length(object$argvals), "colonnes."))
  # }

  # 1. Décomposition des nouvelles données avec la base du modèle
  if (object$method == "bsplines") {
    fd_new <- fda::Data2fd(argvals = object$argvals, y = t(newdata), basisobj = object$basis)
    new_coefs <- t(fd_new$coefs)
  } else if (object$method == "wavelets") {
    # On applique le même filtre et le même niveau qu'à l'entraînement !
    new_coefs <- t(apply(newdata, 1, function(x) {
      wt_x <- wavelets::dwt(x, filter = object$filter, n.levels = object$n_levels, boundary = object$boundary)
      return(wt_x@W[[object$level]])
    }))
  }

  # 2. Projection PLS et LDA
  pls_scores <- predict(object$pls_model, comps = 1:object$ncomp, newdata = new_coefs, type = "scores")

  # 3. Filtrage des composantes constantes
  clean_scores <- data.frame(pls_scores[, !object$const_idx])

  # On force les noms de colonnes à être identiques à ceux du modèle LDA original
  colnames(clean_scores) <- colnames(object$lda_model$means)

  lda_res <- predict(object$lda_model, newdata = clean_scores)

  # 3. Application du seuil
  final_class <- as.character(lda_res$class)
  max_probs <- apply(lda_res$posterior, 1, max)

  if (!is.null(threshold)) {
    final_class[max_probs < threshold] <- "Indétermine"
  }

  return(list(
    class = final_class,
    Probabilite = max_probs,
    lda_res = lda_res
    # check.names = FALSE
  ))
}



#' Resume statistique du modele
#' @param object ...
#' @export
summary_spectra <- function(object) {
  cat("Résumé du Modèle maldiscrim \n")
  cat("Méthode de base :", object$method, "\n")
  if(object$method == "bsplines") cat("Nombre de bases  :", object$nbasis, "\n")
  cat("Composantes PLS  :", object$ncomp, " (", sum(object$const_idx), " supprimées)\n", sep="")
  cat("Nombre de souches :", length(object$labels), "\n")
  cat("Souches identifiées :", paste(object$labels, collapse=","), "\n")
}
