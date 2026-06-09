## code to prepare `ignore_models` dataset goes here

# Charger les donnees du package
data("spectra100", package = "maldiscrim")

# Creer le modele
fpls_model <- fitFPLS_DA(data = spectra100, method = "bsplines", nbasis = 1050)

# L'enregistrer dans sysdata.rda
usethis::use_data(fpls_model, overwrite = TRUE)
