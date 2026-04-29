#' Importation et Pretraitement des replicas biologiques de souches
#'
#' @param data_file Chemin vers le dossier contenant les sous-dossiers de souches.
#' @param method_baseline Methode de retrait de la ligne de base ("SNIP" par defaut).
#' @param halfWindowSize Fenetre de valeurs pour le lissage SavitzkyGolay.
#' @return Une matrice d'intensites (MSP) ou chaque ligne est un replica biologique normalise.
#' @export
process_maldi <- function(data_file, method_baseline = "SNIP", halfWindowSize = 10) { # length intensity à ajouter

  # 1. Liste des dossiers de souches
  souches <- list.dirs(data_file, full.names = FALSE, recursive = FALSE)

  rep_bio <- do.call(rbind, lapply(souches, function(souche){
    file_souche <- file.path(data_file, souche)

    ### Importation via MALDIquantForeign
    spectra <- MALDIquantForeign::import(file_souche, verbose = FALSE)

    ### Pipeline de prétraitement
    #Stabilisation de la variance
    spectra <- MALDIquant::transformIntensity(spectra, method = "sqrt")
    #Smoothing
    spectra <- MALDIquant::smoothIntensity(spectra, method = "SavitzkyGolay", halfWindowSize = halfWindowSize)
    #Suppression de la ligne de base
    spectra <- MALDIquant::removeBaseline(spectra, method = method_baseline, iterations = 100)
    #Etalonnage/normalisation de l'intensite
    spectra <- MALDIquant::calibrateIntensity(spectra, method = "TIC")

    # Groupement par Réplicas Biologiques (via métadonnées). On utilise la logique Nom du réplicas + Mois
    strain <- factor(sapply(spectra, function(x) {
      paste0(MALDIquant::metaData(x)$fullName, "_", lubridate::month(MALDIquant::metaData(x)$acquisitionDate))
    }))
    strain_spectra <- split(spectra, strain)


    # Alignement et calcul du spectre moyen (MSP) pour chaque réplica bio
    msp_list <- lapply(strain_spectra, function(spectralist) {
      aligned <- MALDIquant::alignSpectra(spectralist, halfWindowSize = 20, SNR = 2, tolerance = 0.001, warpingMethod = "lowess")
      return(MALDIquant::averageMassSpectra(aligned, method = "mean"))
    })


    # Chaque colonne est un replica bio (MSP) et on transforme en vecteur numérique (longueur fixé à 20664 par défaut)
    msp_souche <- do.call(cbind, lapply(msp_list, function(msp_obj) {
      intensity <- msp_obj@intensity
      length(intensity) <- 20664
      intensity[is.na(intensity)] <- 0
      return(intensity)
    }))

    # On transpose la matrice msp_souche (lignes = réplicats bio, colonnes = masses)
    msp_souche <- t(msp_souche)

    # On identifie les lignes avec le nom du dossier parent (la souche)
    rownames(msp_souche) <- rep(souche, nrow(msp_souche))

    return(msp_souche)
  }))

  return(rep_bio)
}
