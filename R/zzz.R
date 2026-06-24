.onLoad <- function(libname, pkgname) {

  # Disable reticulate >= 1.41 auto-managed ephemeral environments (uv).
  # Without this, any call that triggers Python initialisation creates a
  # temporary uv env that partially loads TensorFlow, which then conflicts
  # with our own maldiscrim-env when it is later activated.
  Sys.setenv(RETICULATE_USE_MANAGED_VENV = "no")

  # Store virtual environments in a package-specific cache folder that is
  # guaranteed to be outside cloud-synced directories such as OneDrive.
  venvHome <- tools::R_user_dir("maldiscrim", which = "cache")
  if (!dir.exists(venvHome)) dir.create(venvHome, recursive = TRUE)
  Sys.setenv(WORKON_HOME = venvHome)

  # If maldiscrim-env already exists, pre-configure reticulate to use it
  # BEFORE Python is initialised. use_virtualenv(..., required = FALSE)
  # registers the preference without starting a Python session, so it takes
  # priority over the auto-configure fallback. This means the user does not
  # need to call maldiscrim_install_python() at the start of every R session.
  venvPath <- file.path(venvHome, "maldiscrim-env")
  if (dir.exists(venvPath)) {
    reticulate::use_virtualenv(venvPath, required = FALSE)
  }



  # Solve FNN pre-training model paths
#   fnn_env <- tryCatch(
#     get("fnn_model", envir = asNamespace(pkgname)),
#     error = function(e) NULL
#   )
#   if (!is.null(fnn_env)) {
#     fnn_env$modelPath <- system.file(fnn_env$modelPath, package = pkgname)
#     fnn_env$dataPath  <- system.file(fnn_env$dataPath,  package = pkgname)
#   }


  tryCatch({
    env <- asNamespace(pkgname)
    if (exists("fnn_model", envir = env)) {
      obj <- get("fnn_model", envir = env)
      resolve <- function(path) {
        if (!is.null(path) && startsWith(path, ":package:")) {
          relative <- sub("^:package:", "", path)
          system.file(relative, package = pkgname)
        } else {
          path
        }
      }
      obj$modelPath <- resolve(obj$modelPath)
      obj$dataPath  <- resolve(obj$dataPath)
      assign("fnn_model", obj, envir = env)
    }
  }, error = function(e) NULL)


 }

