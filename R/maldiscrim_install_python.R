#' Install the Python Environment Required for the FNN Method
#'
#' @description
#' Creates a dedicated Python virtual environment named \code{"maldiscrim-env"} and installs
#' all Python packages required by the functional neural network (FNN) part of the package.
#'
#' @param pythonVersion Character. Python version used to create the virtual environment. Default is \code{"3.10"}.
#' @param force Logical. If \code{TRUE}, removes and recreates the environment if it already exists. Default is \code{FALSE}.
#' @param verbose Logical. If \code{TRUE}, prints progress messages. Default is \code{TRUE}.
#'
#' @details
#' This function only needs to be run once per machine. It installs \code{tensorflow}, \code{keras}, \code{numpy},
#' \code{pandas}, \code{scikit-learn}, \code{scipy}, \code{umap-learn}, \code{plotly}, \code{matplotlib}, \code{seaborn},
#' \code{plotnine}, \code{opencv-python}, \code{tf-keras-vis} and \code{xplique} into the \code{"maldiscrim-env"}
#' virtual environment, then activates it for the current R session via \code{reticulate::use_virtualenv()}.
#'
#' The environment is stored in a package-specific cache folder rather than the default location used by
#' \code{reticulate}, to avoid failures on machines where that default folder is inside a cloud-synced directory
#' such as OneDrive.
#'
#' @return Invisibly returns \code{TRUE} if the environment was successfully created or already existed and was activated.
#'
#' @seealso \code{\link{fitFNN}}
#'
#' @examples
#' \dontrun{
#' maldiscrim_install_python()
#' }
#'
#' @export
maldiscrim_install_python <- function(pythonVersion = "3.10", force = FALSE, verbose = TRUE) {
  if (!is.character(pythonVersion) || length(pythonVersion) != 1) {
    stop("'pythonVersion' must be a single character string.")
  }
  if (!is.logical(force) || length(force) != 1) {
    stop("'force' must be a single logical value.")
  }
  envName <- "maldiscrim-env"
  .setMaldiscrimVenvHome()
  pythonPackages <- c(
    "tensorflow>=2.12,<2.16", "numpy < 2", "pandas", "scikit-learn", "scipy",
    "umap-learn", "plotly", "matplotlib", "seaborn", "plotnine", "opencv-python", "tf-keras-vis", "xplique"
  )
  envExists <- envName %in% reticulate::virtualenv_list()
  if (envExists && force) {
    if (verbose) message(sprintf("Removing existing environment '%s'...", envName))
    reticulate::virtualenv_remove(envName, confirm = FALSE)
    envExists <- FALSE
  }
  if (!envExists) {
    if (verbose) message(sprintf("Creating virtual environment '%s'...", envName))
    reticulate::virtualenv_create(
      envname  = envName,
      version  = pythonVersion,
      packages = pythonPackages
    )
    if (verbose) message("Environment created and packages installed.")
  } else {
    if (verbose) message(sprintf("Environment '%s' already exists. Use force = TRUE to recreate it.", envName))
  }
  reticulate::use_virtualenv(envName, required = TRUE)
  if (verbose) message(sprintf("Environment '%s' activated for this R session.", envName))
  invisible(TRUE)
}
