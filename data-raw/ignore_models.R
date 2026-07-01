## code to prepare `ignore_models` dataset goes here

# Charger les donnees du package
data("spectra100", package = "maldiscrim")

# Creer le modele
fpls_model <- fitFPLS_DA(data = spectra100, method = "bsplines", nbasis = 1050)

# L'enregistrer dans sysdata.rda
usethis::use_data(fpls_model, overwrite = TRUE)


#### -------------Save pre_training FNN model-------------------------
library(reticulate)
# S'assurer qu'on est dans le bon environnement
reticulate::use_virtualenv("maldiscrim-env", required = TRUE)

tf <- import("tensorflow")

# Charger le modèle actuel
loaded_model <- tf$keras$models$load_model("outputs/fnn/trainedFNN/fnn.keras")


# Créer les dossiers dans inst/ -----------------------------------------
dir.create("inst/models",  recursive = TRUE, showWarnings = FALSE)
dir.create("inst/extdata", recursive = TRUE, showWarnings = FALSE)

# Convertir et sauvegarder en .h5 ---------------------------------------
loaded_model$save("inst/models/fnn_base_model.h5")

# Copier le CSV ---------------------------------------------------------
file.copy(model$dataPath, "inst/extdata/FNN_input.csv")

# Vérifier que les deux fichiers sont bien là ---------------------------
file.exists("inst/models/fnn_base_model.h5")
file.exists("inst/extdata/FNN_input.csv")

# Mettre à jour les chemins et sauvegarder l'objet ----------------------
fnn_model <- model
fnn_model$modelPath  <- ":package:models/fnn_base_model.h5"
fnn_model$dataPath   <- ":package:extdata/FNN_input.csv"
fnn_model$outputPath <- NULL
usethis::use_data(fnn_model, internal = FALSE, overwrite = TRUE)


devtools::install(upgrade = "never")
library(maldiscrim)
data("fnn_model")
cat("modelPath :", fnn_model$modelPath, "\n")
cat("dataPath  :", fnn_model$dataPath,  "\n")
file.exists(fnn_model$modelPath)
file.exists(fnn_model$dataPath)
