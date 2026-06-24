#' Interpretability and attribution functions for FNN models
#'
#' @description
#' This suite of functions provides diagnostic tools to interpret a fitted \code{FNN} model
#' by identifying influential spectral regions (m/z values) using Grad-CAM, Occlusion Sensitivity,
#' first-layer Convolutional Filter Activations, and Activation Maximization.
#'
#' @param object A fitted model object of class inheriting from \code{"FNN"}.
#' @param classIndex Integer. The 1-based index of the target class to compute attributions or maximization for.
#' @param sampleIndex Integer. The index of the specific sample in the dataset to visualize. Default is \code{1}.
#' @param windowSize Integer. The size of the occlusion window in data steps. Default is \code{10}.
#' @param steps Integer. The number of gradient ascent optimization steps for activation maximization. Default is \code{50}.
#' @param palette Either \code{NULL} (automatic) or a user-supplied palette: a character vector of colours or a named \code{RColorBrewer} palette string.
#' @param ... Further arguments passed to or from other internal methods.
#'
#' @details
#' \code{plotGradCamFNN} computes activation heatmaps from the last functional convolutional layer
#' to highlight regions driving class decisions.
#'
#' \code{plotOcclusionFNN} systematically perturbs spectral intervals to monitor performance drops,
#' mapping sensitivity profiles across the m/z domain.
#'
#' \code{plotActivationFNN} extracts and visualizes the transformed feature map outputs generated
#' by each filter within the first functional convolutional layer.
#'
#' \code{plotActivationMaxFNN} uses gradient ascent to synthesize an ideal pseudo-spectrum
#' that maximizes the prediction score for a specific target class.
#'
#' @examples
#' \donttest{
#' # Assuming 'fnn_model' is successfully trained via fitFNN
#'
#' # 1. Plot Grad-CAM attributions
#' plotGradCamFNN(fnn_model, classIndex = 1, sampleIndex = 1)
#'
#' # 2. Plot Occlusion Sensitivity profile
#' plotOcclusionFNN(fnn_model, classIndex = 1, sampleIndex = 1, windowSize = 20)
#'
#' # 3. Plot First-Layer Convolutional Filter Activations
#' plotActivationFNN(fnn_model, sampleIndex = 1)
#'
#' # 4. Plot Synthesized Class Activation Maximization Profile
#' plotActivationMaxFNN(fnn_model, classIndex = 1, steps = 50)
#' }
#'
#' @name plotXplainFNN
NULL

# ══════════════════════════════════════════════════════════════════════════════
# 1. plotGradCamFNN
# ══════════════════════════════════════════════════════════════════════════════

#' @rdname plotXplainFNN
#' @importFrom ggplot2 ggplot aes geom_line geom_rect scale_fill_gradient
#'   theme_minimal labs theme element_blank
#' @export
plotGradCamFNN <- function(object, classIndex, sampleIndex = 1, ...) {

  if (!inherits(object, "FNN")) stop("'object' must be a fitted FNN model.")
  if (missing(classIndex))      stop("'classIndex' must be specified.")

  object <- .resolveFNNPaths(object)

  .checkPythonEnv()

  tf    <- .tf()
  keras <- .keras()

  fnnModel  <- .loadFNNModel(object$modelPath)
  data_list <- .loadFeaturesFromCSV(object)
  features  <- data_list$features

  py_sample_idx <- as.integer(sampleIndex - 1L)
  py_class_idx  <- as.integer(classIndex  - 1L)
  sample_tensor <- tf$expand_dims(features[py_sample_idx], 0L)

  # Locate the last FunctionalConvolution layer
  layers_list     <- fnnModel$layers
  conv_layer_name <- NULL
  for (i in rev(seq_along(layers_list))) {
    if (grepl("FunctionalConvolution", layers_list[[i]]$name)) {
      conv_layer_name <- layers_list[[i]]$name
      break
    }
  }
  if (is.null(conv_layer_name)) {
    stop("No FunctionalConvolution layer found in model configuration.")
  }

  conv_layer <- fnnModel$get_layer(conv_layer_name)

  grad_model <- keras$models$Model(
    inputs  = fnnModel$input,
    outputs = reticulate::tuple(conv_layer$output, fnnModel$output)
  )

  # GradientTape via reticulate::with — replaces the broken %as% syntax
  tape <- tf$GradientTape()
  reticulate::with(tape, {
    res         <- grad_model(sample_tensor)
    conv_outs   <- res[[1]]
    predictions <- res[[2]]
    loss_value  <- predictions[0L, py_class_idx]
  })

  grads        <- tape$gradient(loss_value, conv_outs)
  pooled_grads <- tf$reduce_mean(grads, axis = reticulate::tuple(0L, 1L))

  c_outs  <- reticulate::py_to_r(conv_outs[0L])
  p_grads <- reticulate::py_to_r(pooled_grads)

  # Weighted sum of feature maps — vectorised, no loop
  heatmap <- as.vector(c_outs %*% p_grads)
  heatmap <- pmax(heatmap, 0)
  if (max(heatmap) > 0) heatmap <- heatmap / max(heatmap)

  raw_spectrum <- reticulate::py_to_r(features[py_sample_idx, , 1])
  n_features   <- length(raw_spectrum)
  x_steps      <- seq_len(n_features)

  heatmap_interp <- stats::approx(
    x    = seq(1, n_features, length.out = length(heatmap)),
    y    = heatmap,
    xout = x_steps
  )$y

  df_plot <- data.frame(
    Index      = x_steps,
    Intensity  = raw_spectrum,
    Importance = heatmap_interp
  )

  target_class_name <- object$classNames[classIndex]

  p <- ggplot2::ggplot(df_plot) +
    ggplot2::geom_rect(
      ggplot2::aes(xmin = Index - 0.5, xmax = Index + 0.5,
                   ymin = -Inf, ymax = Inf, fill = Importance),
      alpha = 0.3
    ) +
    ggplot2::geom_line(ggplot2::aes(x = Index, y = Intensity),
                       color = "#377EB8", linewidth = 0.8) +
    ggplot2::scale_fill_gradient(low = "white", high = "#E41A1C") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title    = sprintf("Grad-CAM Spectrum Attribution (Sample %d)",
                         sampleIndex),
      subtitle = sprintf("Target Evaluation Class: %s", target_class_name),
      x        = "Spectral Coordinate Index (m/z)",
      y        = "Intensity Profile",
      fill     = "Activation"
    ) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

  return(p)
}


# plotOcclusionFNN ------------------------------------------------------------

#' @rdname plotXplainFNN
#' @importFrom ggplot2 ggplot aes geom_line geom_rect scale_fill_gradient
#'   theme_minimal labs theme element_blank
#' @export
plotOcclusionFNN <- function(object, classIndex, sampleIndex = 1,
                             windowSize = 10, ...) {

  if (!inherits(object, "FNN")) stop("'object' must be a fitted FNN model.")
  if (missing(classIndex))      stop("'classIndex' must be specified.")

  object <- .resolveFNNPaths(object)

  .checkPythonEnv()

  tf <- .tf()

  fnnModel   <- .loadFNNModel(object$modelPath)
  data_list  <- .loadFeaturesFromCSV(object)
  features   <- data_list$features
  n_features <- data_list$nFeatures

  py_sample_idx <- as.integer(sampleIndex - 1L)
  py_class_idx  <- as.integer(classIndex  - 1L)

  # Extract the sample as an R matrix BEFORE creating any tensor
  # — TF tensors are immutable, R matrices are not
  sample_matrix <- reticulate::py_to_r(features[py_sample_idx])

  baseline_tensor <- tf$expand_dims(
    tf$constant(sample_matrix, dtype = tf$float32), 0L
  )
  baseline_preds <- fnnModel$predict(baseline_tensor, verbose = 0L)
  baseline_prob  <- reticulate::py_to_r(baseline_preds)[1L, py_class_idx + 1L]

  window_starts <- seq(1L, n_features, by = as.integer(windowSize))

  importance_scores <- vapply(window_starts, function(start_pos) {
    end_pos <- min(start_pos + windowSize - 1L, n_features)

    # Mask on the R matrix — safe, no TF item assignment
    masked_matrix <- sample_matrix
    masked_matrix[start_pos:end_pos, 1L] <- 0

    masked_tensor <- tf$expand_dims(
      tf$constant(masked_matrix, dtype = tf$float32), 0L
    )
    masked_preds <- fnnModel$predict(masked_tensor, verbose = 0L)
    masked_prob  <- reticulate::py_to_r(masked_preds)[1L, py_class_idx + 1L]

    max(0, baseline_prob - masked_prob)
  }, numeric(1))

  # Expand window scores back to full feature resolution
  full_importance <- numeric(n_features)
  for (i in seq_along(window_starts)) {
    start_pos <- window_starts[i]
    end_pos   <- min(start_pos + windowSize - 1L, n_features)
    full_importance[start_pos:end_pos] <- importance_scores[i]
  }
  if (max(full_importance) > 0) {
    full_importance <- full_importance / max(full_importance)
  }

  df_plot <- data.frame(
    Index      = seq_len(n_features),
    Intensity  = sample_matrix[, 1L],
    Importance = full_importance
  )

  target_class_name <- object$classNames[classIndex]

  p <- ggplot2::ggplot(df_plot) +
    ggplot2::geom_rect(
      ggplot2::aes(xmin = Index - 0.5, xmax = Index + 0.5,
                   ymin = -Inf, ymax = Inf, fill = Importance),
      alpha = 0.3
    ) +
    ggplot2::geom_line(ggplot2::aes(x = Index, y = Intensity),
                       color = "#377EB8", linewidth = 0.8) +
    ggplot2::scale_fill_gradient(low = "white", high = "#FF7F00") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title    = sprintf("Occlusion Sensitivity Profile (Sample %d)",
                         sampleIndex),
      subtitle = sprintf("Target Evaluation Class: %s (Window Size: %d)",
                         target_class_name, windowSize),
      x        = "Spectral Coordinate Index (m/z)",
      y        = "Intensity Profile",
      fill     = "Sensitivity"
    ) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

  return(p)
}


# plotActivationFNN -----------------------------------------------------------

#' @rdname plotXplainFNN
#' @importFrom ggplot2 ggplot aes geom_line facet_wrap theme_minimal labs
#' @export
plotActivationFNN <- function(object, sampleIndex = 1, ...) {

  if (!inherits(object, "FNN")) stop("'object' must be a fitted FNN model.")

  object <- .resolveFNNPaths(object)

  .checkPythonEnv()

  # Bug fix: tf was missing — keras alone was not sufficient
  tf    <- .tf()
  keras <- .keras()

  fnnModel  <- .loadFNNModel(object$modelPath)
  data_list <- .loadFeaturesFromCSV(object)
  features  <- data_list$features

  py_sample_idx <- as.integer(sampleIndex - 1L)
  sample_tensor <- tf$expand_dims(features[py_sample_idx], 0L)

  # Locate the first FunctionalConvolution layer
  first_conv_name <- NULL
  for (layer in fnnModel$layers) {
    if (grepl("FunctionalConvolution", layer$name)) {
      first_conv_name <- layer$name
      break
    }
  }
  if (is.null(first_conv_name)) {
    stop("No functional convolutional layer found in model structure.")
  }

  sub_model <- keras$models$Model(
    inputs  = fnnModel$input,
    outputs = fnnModel$get_layer(first_conv_name)$output
  )
  activations_tensor <- sub_model$predict(sample_tensor, verbose = 0L)
  activations_matrix <- reticulate::py_to_r(activations_tensor[0L])

  n_features <- nrow(activations_matrix)
  n_filters  <- ncol(activations_matrix)

  # Vectorised construction — no loop, no do.call
  df_plot <- data.frame(
    Index      = rep(seq_len(n_features), times = n_filters),
    Activation = as.vector(activations_matrix),
    Filter     = rep(sprintf("Filter %02d", seq_len(n_filters)),
                     each = n_features)
  )

  p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = Index, y = Activation)) +
    ggplot2::geom_line(color = "#4DAF4A", linewidth = 0.7) +
    ggplot2::facet_wrap(~Filter, scales = "free_y") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title    = "First-Layer Convolutional Filter Activations",
      subtitle = sprintf("Evaluated on Data Sample Index: %d", sampleIndex),
      x        = "Feature Coordinate Index",
      y        = "Transformation Output Amplitude"
    )

  return(p)
}


# plotActivationMaxFNN --------------------------------------------------------

#' @rdname plotXplainFNN
#' @importFrom ggplot2 ggplot aes geom_line theme_minimal labs
#' @export
plotActivationMaxFNN <- function(object, classIndex, steps = 50, ...) {

  if (!inherits(object, "FNN")) stop("'object' must be a fitted FNN model.")
  if (missing(classIndex))      stop("'classIndex' must be specified.")

  object <- .resolveFNNPaths(object)

  .checkPythonEnv()

  tf <- .tf()

  fnnModel     <- .loadFNNModel(object$modelPath)
  input_shape  <- fnnModel$input$shape
  n_features   <- as.integer(input_shape[1])
  py_class_idx <- as.integer(classIndex - 1L)

  # Initialise pseudo-spectrum from uniform random noise bounded to [0, 0.1]
  input_data <- tf$Variable(
    tf$random$uniform(
      shape   = reticulate::tuple(1L, n_features, 1L),
      minval  = 0.0,
      maxval  = 0.1
    )
  )

  # Gradient ascent — reticulate::with replaces broken %as% syntax
  for (step in seq_len(steps)) {
    tape <- tf$GradientTape()
    reticulate::with(tape, {
      tape$watch(input_data)
      predictions <- fnnModel(input_data)
      loss_value  <- predictions[0L, py_class_idx]
    })

    grads <- tape$gradient(loss_value, input_data)
    grads <- grads / (tf$math$reduce_std(grads) + 1e-5)
    input_data$assign_add(grads * 0.01)
    input_data$assign(tf$clip_by_value(input_data, 0.0, 1.0))
  }

  optimized_spectrum <- reticulate::py_to_r(input_data)[1L, , 1L]

  df_plot <- data.frame(
    Index     = seq_along(optimized_spectrum),
    Intensity = optimized_spectrum
  )

  target_class_name <- object$classNames[classIndex]

  p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = Index, y = Intensity)) +
    ggplot2::geom_line(color = "#984EA3", linewidth = 0.9) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title    = "Class Activation Maximization Profile",
      subtitle = sprintf(
        "Synthesized Idealized Spectrum for Target Class: %s",
        target_class_name
      ),
      x = "Spectral Coordinate Index (m/z)",
      y = "Optimized Feature Amplitude"
    )

  return(p)
}
