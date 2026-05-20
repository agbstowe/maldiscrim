#' Fit Functionnal Partiel Least Square and Discriminant Analysis (fPLS_DA) Model
#'
#' @description
#' `fit.fPLS_DA` performs a supervised classification of MALDI-TOF mass spectra by combining functional data analysis (FDA)
#' dimensionality reduction techniques with Partial Least Squares Regression (PLSR) and Linear Discriminant Analysis (LDA).
#'
#'
#' @param data A numeric matrix of processed Mean Mass Spectra (MSP) where rows
#' represent biological replicates (with strain names as rownames) and columns represent spectral intensities.
#' @param method Character. The mathematical basis used for spectral decomposition: either `"bsplines"`
#' (B-spline basis functions) or `"wavelets"` (Discrete Wavelet Transform).
#' @param nbasis Integer. Number of basis functions used when `method = "bsplines"`. Default is `1050`.
#' @param ncomp Integer. Number of PLS components to retain. If `NULL` (default), the optimal number of
#' components is automatically determined using the one-sigma heuristic based on Cross-Validation RMSEP.
#' @param argvals A numeric vector specifying the evaluation grid (m/z indices) for the B-spline projection.
#' @param rangeval A numeric vector of length 2 defining the domain range for the B-spline basis.
#' @param n_folds Integer. Number of folds for the internal cross-validation performed during PLS fitting. Default is `5`.
#' @param filter Character. The wavelet filter type passed to \code{\link[wavelets]{dwt}} when `method = "wavelets"`. Default is `"la14"`.
#' @param boundary Character. The boundary handling method passed to \code{\link[wavelets]{dwt}}. Default is `"periodic"`.
#' @param n.levels Integer. Total number of decomposition levels for the wavelet transform.
#' If `NULL`, it is automatically set to the maximum possible level.
#' @param level Integer. The specific wavelet decomposition resolution level from which
#' coefficients are extracted for downstream classification. Default is `3`.
#' @param precomputed Logical. If `TRUE`, the function bypasses all heavy computations and returns a built-in
#' pre-trained model object. Primarily used to speed up examples and unit testing. Default is `FALSE`.
#'
#'
#' @details
#' The pipeline executes three sequential phases to achieve classification of functional spectral data:
#' \enumerate{
#'   \item **Functional Decomposition:** Depending on `method`, the high-dimensional spectral matrix is projected into either:
#'     \itemize{
#'       \item *B-Splines:* Continuous functional representations are built using a B-spline basis matrix via `fda::Data2fd`.
#'       \item *Wavelets:* Multi-resolution analysis is applied via the Discrete Wavelet Transform (`wavelets::dwt`). A built-in guardrail automatically adjusts `level` if the remaining coefficients are insufficient for modeling.
#'     }
#'   \item **PLS Dimension Reduction:** A Kernel Partial Least Squares Regression (`pls::plsr`) is fitted using dummy-coded groups as responses.
#'   \item **Classification:** Low-variance components are filtered out, and a Linear Discriminant Analysis (`MASS::lda`) is fitted on the remaining PLS scores to generate the final decision boundaries.
#' }
#'
#' @return An object of class `"fPLS_DA"`. This object is a named list containing:
#' \item{pls_model}{The fitted PLS model object returned by `pls::plsr`.}
#' \item{lda_model}{The fitted LDA model object returned by `MASS::lda`.}
#' \item{method}{A character string indicating the functional decomposition method used.}
#' \item{basis}{The B-spline basis object (if `method = "bsplines"`), otherwise `NULL`.}
#' \item{wavelets}{A list containing raw DWT outputs for each spectrum (if `method = "wavelets"`), otherwise `NULL`.}
#' \item{ncompOpt}{The optimized or user-defined number of PLS components retained.}
#' \item{const_idx}{A logical vector indicating which PLS components were excluded due to near-zero variance.}
#' \item{levels}{The structural names/categories of the groups (strains).}
#' \item{labels}{The original factor array of group memberships assigned to each row.}
#' \item{argvals}{The grid evaluations used for B-splines projection.}
#' \item{filter}{The name of the wavelet filter used.}
#' \item{n_levels}{The maximum wavelet decomposition level reached.}
#' \item{level}{The specific wavelet resolution level extracted.}
#' \item{nbasis}{The total number of B-spline basis functions used.}
#' \item{boundary}{The wavelet boundary method applied.}
#'
#' @seealso \code{\link[pls]{plsr}}, \code{\link[MASS]{lda}}, \code{\link[fda]{Data2fd}}, \code{\link[wavelets]{dwt}}
#'
#' @examples
#' # Model loading using the precomputed cache
#' model_res <- fit.fPLS_DA(precomputed = TRUE)
#' # print(model_res)
#'
#' \donttest{
#' # Real execution example using the package built-in data 'spectra100'
#' data("spectra100", package = "maldiscrim")
#'
#' # Fit using B-Splines functional basis
#' bspline_model <- fit.fPLS_DA(
#'   data = spectra100, method = "bsplines", nbasis = 1050, n_folds = 5
#'   )
#'
#' # Fit using Wavelets multi-resolution decomposition
#' # wavelet_model <- fit.fPLS_DA(
#' # data = spectra100, method = "wavelets", level = 3, n_folds = 5)
#' }
#'
#' @importFrom stats sd aggregate
#' @export
fit.fPLS_DA <- function(data, method = c("bsplines", "wavelets"), nbasis = 1050, ncomp = NULL, argvals = 2001:22664,
                          rangeval = c(1999,22664), n_folds = 5,
                          filter = "la14", boundary = "periodic", n.levels = NULL, level = 3, precomputed = FALSE) {

  # Interception pour charger le resultat pre-calcule instantanement (Gain de temps exemples/tests)
  if (precomputed) {
    message("Loading precomputed fPLS_DA model from cache")
    # S'assure que l'objet cache existe dans le package
    if (exists("fpls_model", envir = asNamespace("maldiscrim"))) {
      return(get("fpls_model", envir = asNamespace("maldiscrim")))
    } else {
      stop("Cached model 'fpls_model' not found in the package data.")
    }
  }

  method <- match.arg(method)
  groups <- as.factor(rownames(data))
  # if (is.null(argvals)){
  #   argvals <- 2001:22664  ## Grille de projection de tous les spectres. Tous les spectres doivent être calibrés sur cet intervalle
  # } else {
  #   argvals = argvals
  # }

  res_fdaDecomposition <- fdaDecomposition(data,method,nbasis,rangeval,argvals)

  # 1. Décomposition selon la méthode choisie
  if (method == "bsplines") {
    message("Utilisation de la decomposition fonctionnelle en B-splines.")
    basis_obj <- fda::create.bspline.basis(rangeval = rangeval, nbasis = nbasis)
    fd_obj <- fda::Data2fd(argvals = argvals, y = t(data), basisobj = basis_obj)
    coef_matrix <- t(fd_obj$coefs)
    }
  else if (method == "wavelets") {
    message("Utilisation de la decomposition fonctionnelle en ondelettes.")

    max_levels <- floor(log2(ncol(data)))

    if (is.null(n.levels)) {
      # On prend la puissance de 2 la plus proche de la longueur du signal
      # log2(20664) donne environ 14.3, donc floor() nous donne 14
      n.levels <- max_levels
    } else if (n.levels > max_levels) {
      warning(paste("n.levels trop grand. Regle automatiquement sur", max_levels))
      n.levels <- max_levels
    }

    # 3. Validation de la pertinence du 'level' (Garde-fou scientifique)
    # On s'assure qu'au niveau choisi, il reste assez de coefficients (ici > 20 pour éviter que la PLS ne plante)
    # Nombre theorique de coeffs ~= Longueur / 2^level
    n_coef_theorique <- floor(ncol(data) / (2^level))

    if (n_coef_theorique < 20) {
      new_level <- max(1, floor(log2(ncol(data) / 20)))
      warning(paste("Le level choisi ", level, " est trop eleve pour satisfaire la PLS. Ajustement automatique a : ", new_level))
      level <- new_level
    }

    if (level > n.levels) {
      warning(paste("Le level choisi ", level, " superieur au nombre total de niveaux. Ajustement a ", n.levels))
      level <- n.levels
    }

    # On applique la DWT sur chaque spectre (chaque ligne de data)
    # On utilise 'boundary="periodic"' pour eviter les effets de bord au début/fin du spectre
    # On utilise apply pour creer une liste contenant les objets dwt
    wt_list <- apply(data, 1, function(x) {
      wavelets::dwt(x, filter = filter, n.levels = n.levels, boundary = boundary)
    })

    # Extraction des coefficients pour la PLS
    coef_matrix <- t(sapply(wt_list, function(wt) wt@W[[level]]))
    wt_output <- wt_list # On stocke la liste complète
  }

  # 2. PLS sur ncomp composantes et Gestion de la variance nulle
  y_dummy <- fastDummies::dummy_cols(groups, remove_selected_columns = TRUE)

  max_search <- min(nrow(data) - 1, 20) # On cherche jusqu'a 20 composantes max pour la pls pour eviter un eventuelle sur apprentissage

  pls_mod <- pls::plsr(as.matrix(y_dummy) ~ coef_matrix, ncomp = max_search, method = "kernelpls",
                       validation = "CV", segments = n_folds)

  # Optimisation automatique du nombre de composante
  if (!is.null(ncomp)) {
    # L'utilisateur a fourni une valeur pour ncomp et on borne la valeur par max_search pour éviter les plantages.
    ncomp <- min(ncomp, max_search)
    message("Nombre de composantes optimise automatiquement a ", ncomp, " pour la PLS")

  } else {
    # Calcul automatique du ncomp optimal via la methode onesigma de la fonction selectNcomp
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
    message("Nombre de composantes optimise automatiquement a ", ncomp, " pour la PLS")
  }

  # extraction
  scores <- pls_mod$scores[, 1:ncomp]
  
  # Calcul du filtrage de variance sur les scores
  group_means <- aggregate(scores, list(groups), mean)[,-1]
  f1 <- apply(scores - as.matrix(group_means[groups, ]), 2, sd)
  const_idx <- f1 < 0.0001
  scores_clean <- scores[, !const_idx]

  # 3. Modele LDA
  lda_mod <- MASS::lda(x = scores_clean, grouping = groups)

  # 4. Construction de l'objet S3
  model_obj <- list(
    pls_model = pls_mod,
    lda_model = lda_mod,
    method = method,
    basis = if(method == "bsplines") basis_obj else NULL,
    wavelets = if(method == "wavelets") wt_output else NULL,
    ncompOpt = ncomp,
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

#' Predict Method for fPLS_DA Model Fits
#'
#' @description
#' Obtains predictions (classified strains and posterior probabilities) from a fitted \code{fPLS_DA} model object.
#'
#' @method predict fPLS_DA
#'
#' @param object A fitted model object of class \code{"fPLS_DA"}.
#' @param newdata An optional matrix or data frame of new MALDI-TOF spectral intensities to predict.
#' If omitted (\code{NULL}), the function automatically retrieves the coefficients from the training data and returns the fitted predictions.
#' @param threshold Numeric. A classification probability threshold (between 0 and 1).
#' If the maximum posterior probability for a sample is below this threshold, its predicted class is labeled as \code{"Doubty"}.
#' Default is \code{NULL} (no threshold applied).
#' @param ... Further arguments passed to or from other methods.
#' @importFrom stats predict
#'
#' @details
#' When new spectral data is supplied via \code{newdata}, the function projects it into the exact same functional domain
#' space established during training (using either the original B-spline basis parameters or the specific wavelet filter and resolution level).
#'
#' The extracted functional coefficients are then projected onto the optimal Partial Least Squares (PLS) components.
#' Finally, after removing near-zero variance components, the Linear Discriminant Analysis boundary profiles
#' evaluate the scores to output the most probable strain class along with its associated posterior probability.
#'
#' @return A named list containing three components:
#' \item{class}{A character vector of predicted class/strain assignments for each observation.}
#' \item{Probability}{A numeric vector indicating the maximum posterior probability (confidence score) associated with each prediction.}
#' \item{lda_res}{The raw list output from the underlying \code{\link[MASS]{predict.lda}} step.}
#'
#' @seealso \code{\link{fit.fPLS_DA}}, \code{\link[MASS]{predict.lda}}, \code{\link[pls]{predict.mvr}}
#'
#'
#' @examples
#' # Load the pre-trained example model from the package cache
#' model_res <- fit.fPLS_DA(precomputed = TRUE)
#'
#' # Predict on the training data by omitting 'newdata' (Fitted values)
#' predict_res <- predict(model_res)
#' table(predict_res$class)
#'
#' \donttest{
#' # Example with 'newdata' using the built-in dataset
#' data("spectra100", package = "maldiscrim")
#'
#' # Let's simulate new data by taking a subset of spectra100
#' new_spectra <- spectra100[1:5, ]
#'
#' # Predict
#' predict_res <- predict(model_res, newdata = new_spectra )
#' print(predict_res$class)
#' print(predict_res$Probabilities)
#' }
#'
#' @importFrom stats predict
#' @importFrom MASS lda
#' @import pls
#' @export
predict.fPLS_DA <- function(object, newdata = NULL, threshold = NULL, ...) {

  if (is.null(newdata)) {
    new_coefs <- object$pls_model$model[[2]]

  } else {
    # VALIDATION DEFENSIF
    if (!is.matrix(newdata) && !is.data.frame(newdata)) {
      stop("newdata doit etre une matrice ou un data.frame.")
    }
    # # Vérification cruciale : la taille du spectre. On compare avec la longueur de argvals
    # if (ncol(newdata) != length(object$argvals)) {
    #   stop(paste("Erreur : newdata possede", ncol(newdata),
    #              "colonnes, mais le modèle a ete entraine avec", length(object$argvals), "colonnes."))
    # }

    # 1. Decomposition des nouvelles donnees avec la base du modele
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
  }

  # 2. Projection PLS et LDA
  pls_scores <- predict(object$pls_model, comps = 1:object$ncompOpt, newdata = as.matrix(new_coefs), type = "scores")

  # 3. Filtrage des composantes constantes
  clean_scores <- data.frame(pls_scores[, !object$const_idx])

  # On force les noms de colonnes a etre identiques a ceux du modele LDA original
  colnames(clean_scores) <- colnames(object$lda_model$means)

  lda_res <- predict(object$lda_model, newdata = clean_scores)

  # 3. Application du seuil
  final_class <- as.character(lda_res$class)
  max_probs <- apply(lda_res$posterior, 1, max)

  if (!is.null(threshold)) {
    final_class[max_probs < threshold] <- "Doubty"
  }

  return(list(
    class = final_class,
    Probabilite = max_probs,
    lda_res = lda_res
    # check.names = FALSE
  ))
}



#' Resume statistique du modele
#'
#' @description
#' Description
#'
#' @param object ...
#'
#' @details
#' Additional details...
#'
#' @return function outputs
#'
#' @examples
#' # example code
#'
#'
#' @export
summary_spectra <- function(object) {
  cat("Resume du Modele maldiscrim \n")
  cat("Methode de base :", object$method, "\n")
  if(object$method == "bsplines") cat("Nombre de bases  :", object$nbasis, "\n")
  cat("Composantes PLS  :", object$ncomp, " (", sum(object$const_idx), " supprimees)\n", sep="")
  cat("Nombre de souches :", length(object$labels), "\n")
  cat("Souches identifiees :", paste(object$labels, collapse=","), "\n")
}
