#' Diagnostic visualization functions for FNN outputs
#'
#' @description
#' This group of functions provides a comprehensive suite of tools to visually explore and diagnose the training
#' and latent space projections of a fitted \code{FNN} model.
#'
#' @param object A fitted model object of class inheriting from \code{"FNN"}.
#' @param on Character. The network space to project and plot: either \code{"logits"} or \code{"softmax"}. Default is \code{"logits"}.
#' @param palette Either \code{NULL} (automatic) or a user-supplied palette: a character vector of colours or a named \code{RColorBrewer} palette string.
#' @param seed Integer. Random seed for UMAP and MDS reproducibility. Default is \code{42}.
#' @param ... Further arguments passed to or from other internal methods.
#'
#' @details
#' The baseline performance curves (\code{plotTrainingFNN}) monitor the model convergence and identify potential overfitting.
#' The classification diagnostic layer (\code{plotConfusionFNN}) exposes cross-class ambiguities on the validation set.
#'
#' Manifold embedding plots (\code{plotUMAPFNN} and \code{plotMDSFNN}) evaluate cluster separation across different latent levels
#' (logits vs. softmax activations), adding statistical confidence ellipses to enclose distinct biological strains.
#'
#' @examples
#' \donttest{
#' # Assuming 'fnn_model' is successfully trained via fitFNN
#'
#' # 1. Plot Training History (Loss and Accuracy)
#' plotTrainingFNN(fnn_model)
#'
#' # 2. Plot Confusion Heatmap matrix
#' plotConfusionFNN(fnn_model)
#'
#' # 3. Plot UMAP Projection on Layer Logits
#' plotUMAPFNN(fnn_model, on = "logits")
#'
#' # 4. Plot MDS Projection on Softmax Probabilities
#' plotMDSFNN(fnn_model, on = "softmax", seed = 123)
#' }
#'
#' @name plotDiagnosticsFNN
NULL

# ══════════════════════════════════════════════════════════════════════════════
# 1. plotTrainingFNN
# ══════════════════════════════════════════════════════════════════════════════

#' @rdname plotDiagnosticsFNN
#' @importFrom ggplot2 ggplot aes geom_line geom_point facet_wrap
#'   scale_color_manual theme_minimal labs theme
#' @export
plotTrainingFNN <- function(object, ...) {

  if (!inherits(object, "FNN")) stop("'object' must be a fitted FNN model.")

  history <- object$history
  if (is.null(history)) stop("No training history found in the FNN object.")

  epochs <- seq_along(history$loss)

  df <- rbind(
    data.frame(Epoch = epochs, Value = unlist(history$loss),
               Type = "Train",      Metric = "Loss"),
    data.frame(Epoch = epochs, Value = unlist(history$val_loss),
               Type = "Validation", Metric = "Loss"),
    data.frame(Epoch = epochs, Value = unlist(history$accuracy),
               Type = "Train",      Metric = "Accuracy"),
    data.frame(Epoch = epochs, Value = unlist(history$val_accuracy),
               Type = "Validation", Metric = "Accuracy")
  )

  p <- ggplot2::ggplot(df, ggplot2::aes(x = Epoch, y = Value, color = Type)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_point(size = 1.5) +
    ggplot2::facet_wrap(~Metric, scales = "free_y") +
    ggplot2::scale_color_manual(
      values = c("Train" = "#2196F3", "Validation" = "#F44336")
    ) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "FNN Model Training History",
      x     = "Epochs",
      y     = "Value",
      color = "Dataset"
    ) +
    ggplot2::theme(legend.position = "bottom")

  return(p)
}


# plotConfusionFNN ------------------------------------------------------------

#' @rdname plotDiagnosticsFNN
#' @importFrom ggplot2 ggplot aes geom_tile geom_text scale_fill_gradient
#'   theme_minimal labs
#' @export
plotConfusionFNN <- function(object, ...) {

  if (!inherits(object, "FNN")) stop("'object' must be a fitted FNN model.")

  pred        <- predict.FNN(object, newdata = NULL)
  true_labels <- as.character(object$labels)
  # confusion   <- table(Predicted = pred$class, Actual = true_labels)
  confusion   <- caret::confusionMatrix(as.factor(pred$class), as.factor(true_labels))
  df_conf     <- as.data.frame(as.table(confusion))

  p <- ggplot2::ggplot(df_conf,
                       ggplot2::aes(x = Predicted, y = Actual, fill = Freq)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::geom_text(ggplot2::aes(label = Freq), color = "black", size = 4) +
    ggplot2::scale_fill_gradient(low = "#E3F2FD", high = "#1565C0") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "FNN Confusion Matrix",
      x     = "Predicted Class",
      y     = "True Class",
      fill  = "Observations"
    )

  return(p)
}


# plotUMAPFNN -----------------------------------------------------------------

#' @rdname plotDiagnosticsFNN
#' @importFrom ggplot2 ggplot aes geom_point stat_ellipse scale_color_manual
#'   theme_minimal labs
#' @export
plotUMAPFNN <- function(object, on = c("logits", "softmax"),
                        palette = NULL, seed = 42, ...) {

  if (!inherits(object, "FNN")) stop("'object' must be a fitted FNN model.")
  on <- match.arg(on)

  .checkPythonEnv()

  features_space <- .extractNetworkSpace(object, on = on)

  umap_module <- reticulate::import("umap", convert = FALSE)
  reducer     <- umap_module$UMAP(random_state = as.integer(seed))
  umap_coords <- reticulate::py_to_r(reducer$fit_transform(features_space))

  df_proj <- data.frame(
    UMAP1 = umap_coords[, 1],
    UMAP2 = umap_coords[, 2],
    label = as.factor(object$labels)
  )

  n_groups <- nlevels(df_proj$label)
  pal      <- .getPalette(n_groups, palette)

  p <- ggplot2::ggplot(df_proj,
                       ggplot2::aes(x = UMAP1, y = UMAP2, color = label)) +
    ggplot2::geom_point(alpha = 0.8, size = 2.5) +
    ggplot2::stat_ellipse(type = "t", level = 0.95, linetype = 2,
                          alpha = 0.5) +
    ggplot2::scale_color_manual(values = pal) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = sprintf("UMAP Projection on FNN %s", toupper(on)),
      x     = "UMAP Dimension 1",
      y     = "UMAP Dimension 2",
      color = "Strains"
    )

  return(p)
}


# plotMDSFNN ------------------------------------------------------------------

#' #' @rdname plotDiagnosticsFNN
#' #' @importFrom ggplot2 ggplot aes geom_point stat_ellipse scale_color_manual
#' #'   theme_minimal labs
#' #' @export
#' plotMDSFNN <- function(object, on = c("logits", "softmax"),
#'                        palette = NULL, seed = 42, ...) {
#'
#'   if (!inherits(object, "FNN")) stop("'object' must be a fitted FNN model.")
#'   on <- match.arg(on)
#'
#'   .checkPythonEnv()
#'
#'   features_space   <- .extractNetworkSpace(object, on = on)
#'
#'   sklearn_manifold <- reticulate::import("sklearn.manifold", convert = FALSE)
#'   reducer          <- sklearn_manifold$MDS(random_state = as.integer(seed))
#'   mds_coords       <- reticulate::py_to_r(reducer$fit_transform(features_space))
#'
#'   df_proj <- data.frame(
#'     MDS1  = mds_coords[, 1],
#'     MDS2  = mds_coords[, 2],
#'     label = as.factor(object$labels)
#'   )
#'
#'   n_groups <- nlevels(df_proj$label)
#'   pal      <- .getPalette(n_groups, palette)
#'
#'   p <- ggplot2::ggplot(df_proj,
#'                        ggplot2::aes(x = MDS1, y = MDS2, color = label)) +
#'     ggplot2::geom_point(alpha = 0.8, size = 2.5) +
#'     ggplot2::stat_ellipse(type = "t", level = 0.95, linetype = 2,
#'                           alpha = 0.5) +
#'     ggplot2::scale_color_manual(values = pal) +
#'     ggplot2::theme_minimal() +
#'     ggplot2::labs(
#'       title = sprintf("MDS Projection on FNN %s", toupper(on)),
#'       x     = "MDS Dimension 1",
#'       y     = "MDS Dimension 2",
#'       color = "Strains"
#'     )
#'
#'   return(p)
#' }
