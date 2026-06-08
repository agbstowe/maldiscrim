#' Visualisation functions for maldiscrim outputs
#'
#' @description
#' This group of functions provides a comprehensive suite of tools to visually explore, diagnose, and interpret the results of a fitted \code{FPLS_DA} model.
#' It includes classical 2D/3D score projections, UMAP manifold embeddings, variance scree plots, and spectral loading weight profiles.
#'
#' @param object A fitted model object of class inheriting from \code{"FPLS_DA"}.
#' @param comp A numeric vector specifying the projection axes or component to display.
#' @param show_proba Logical. If \code{TRUE}, sample points are colored by group, and their transparency is dynamically scaled based on their
#' maximum LDA posterior probability tier: 1.0 (\eqn{\ge 0.8}), 0.6 (\eqn{\ge 0.5}), or 0.2 (\eqn{< 0.5}).
#' @param palette Either \code{NULL} (automatic) or a user-supplied palette: a character vector of colours or a named \code{RColorBrewer} palette
#'   string. When \code{NULL}, \code{"Set1"} is used for up to 8 groups and \code{pals} for larger numbers of groups.
#' @param ... Further arguments passed to or from other internal methods.
#'
#' @details
#' The score-based plotting functions (\code{plotSpectra2D}, \code{plotSpectra3D},
#' and \code{plotSpectraUMAP}) assist in assessing cluster separation in different latent spaces.
#' When \code{show_proba = TRUE}, they act as a diagnostic layer to visually isolate samples near classification boundaries.
#'
#' Diagnostic plots (\code{plotPLSvar} and \code{plotPLSweights}) help identify
#' the optimal number of components and track down specific mass-to-charge (\eqn{m/z}) channels.
#'
#' @examples
#' \dontrun{
#' # Assuming 'fit_model' is generated
#' # 1. 2D Score Projection
#' plotSpectra2D(fit_model, comp = c(1, 2), show_proba = TRUE)
#'
#' # 2. Interactive 3D Plot
#' plotSpectra3D(fit_model, comp = c(1, 2, 3))
#'
#' # 3. Variance Scree Plot
#' plotPLSvar(fit_model)
#'
#' # 4. Loading Profiles
#' plotPLSweights(fit_model, comp = 1)
#'
#' # 5. Non-linear UMAP Mapping
#' plotSpectraUMAP(fit_model, n_neighbors = 10, min_dist = 0.05)
#' }
#'
#' @importFrom stats predict
#' @name plot_functions
NULL

#' @rdname plot_functions
#' @section 3D Score Projection:
#' Generates an interactive tridimensional scatter plot of the PLS latent scores using \code{plotly}.
#' Axis titles automatically integrate the percentage of explained variance for each selected component.
#'
#' @export
plotSpectra3D <- function(object, comp = c(1, 2, 3), show_proba = FALSE , palette = NULL, ...) {
  scores_df <- as.data.frame(object$pls_model$scores[, comp])
  colnames(scores_df) <- c("X", "Y", "Z")
  scores_df$strain <- as.factor(object$labels)
  pct <- round(object$pls_model$Xvar / object$pls_model$Xtotvar * 100, 2)

  lda_alpha <- .computeLDAalpha(object)
  scores_df$Alpha <- if (show_proba) lda_alpha else rep(0.8, nrow(scores_df))

  n_groups <- nlevels(scores_df$strain)
  pal      <- .getPalette(n_groups, palette)

  get_lims <- function(vec) {
    r <- range(vec)
    d <- diff(r)
    c(r[1] - 0.1 * d, r[2] + 0.1 * d)
  }
  xlims <- get_lims(scores_df$X)
  ylims <- get_lims(scores_df$Y)
  zlims <- get_lims(scores_df$Z)

  p <- plotly::plot_ly(scores_df, x = ~X, y = ~Y, z = ~Z, color = ~strain,
                  type = 'scatter3d', mode = 'markers', colors = pal,
                  marker = list(size = 3, opacity = ~Alpha,...))
  p <- plotly::layout(p,scene = list(xaxis = list(title = paste0("Comp ", comp[1]," (", pct[comp[1]], "%)") , range = xlims),
                                yaxis = list(title = paste0("Comp ", comp[2]," (", pct[comp[2]], "%)") , range = ylims),
                                zaxis = list(title = paste0("Comp ", comp[3]," (", pct[comp[3]], "%)") , range = zlims)
                                ))

  return(p)
}


#' @rdname plot_functions
#' @section 2D Score Projection:
#' Generates a static two-dimensional scatter plot using \code{ggplot2}.
#'
#' @export
plotSpectra2D <- function(object, comp = c(1, 2), show_proba = FALSE, palette = NULL, ...) {
  # Construction du dataframe de base
  df <- data.frame(object$pls_model$scores[, comp])
  colnames(df) <- c("x", "y")
  df$strain <- as.factor(object$labels)
  pct <- round(object$pls_model$Xvar / object$pls_model$Xtotvar * 100, 2)

  lda_alpha <- .computeLDAalpha(object)
  df$Alpha  <- if (show_proba) lda_alpha else rep(0.8, nrow(df))

  n_groups <- nlevels(df$strain)
  pal      <- .getPalette(n_groups, palette)


  range_x <- range(df$x)
  range_y <- range(df$y)
  xlim <- c(range_x[1] - 0.1 * diff(range_x), range_x[2] + 0.1 * diff(range_x))
  ylim <- c(range_y[1] - 0.1 * diff(range_y ), range_y[2] + 0.1 * diff(range_y))

  p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, color = strain)) +
    ggplot2::geom_point(ggplot2::aes(alpha = Alpha), size = 2) +
    ggplot2::scale_alpha_identity() +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = paste0("Comp ", comp[1], " (", pct[comp[1]], "%)"), y = paste0("Comp ", comp[2], " (", pct[comp[2]], "%)")) +
    ggplot2::scale_color_manual(values = pal) +
    ggplot2::coord_cartesian(xlim = xlim, ylim = ylim)

  return(p)
}



#' @rdname plot_functions
#' @section Explained Variance Scree Plot:
#' Computes and displays a bar chart of the individual variance percentage explained by each functional PLS component.
#'
#' @export
plotPLSvar <- function(object) {
  var_exp <- object$pls_model$Xvar / object$pls_model$Xtotvar * 100
  df <- data.frame(Comp = seq_along(var_exp), Var = var_exp)

  ggplot2::ggplot(df[1:object$ncompOpt,], ggplot2::aes(x = Comp, y = Var)) +
    ggplot2::geom_bar(stat = "identity", fill = "steelblue") +
    ggplot2::labs(title = "Variance explained per PLS component", x = "Component", y = "% Variance")+
    ggplot2::theme_minimal()
}


#' @rdname plot_functions
#' @section Spectral Loadings Profile:
#' Plots the loading weight profiles across variables (mass channels).
#' Useful for identifying specific \eqn{m/z} markers that highly influence class segregation.
#'
#' @export
plotPLSweights <- function(object, comp = 1) {
  loading_values <- object$pls_model$loadings[, comp]
  df <- data.frame(m_z = seq_along(loading_values), Weight = loading_values)

  ggplot2::ggplot(df, ggplot2::aes(x = m_z, y = Weight)) +
    ggplot2::geom_line() +
    ggplot2::labs(title = paste("Loading weights — Component", comp), x = "m/z index", y = "Weight") +
    ggplot2::theme_minimal()
}


#' @rdname plot_functions
#' @section UMAP Non-linear Embedding:
#' Projects the full matrix of high-dimensional PLS scores into a non-linear 2D layout via Uniform Manifold Approximation and Projection (UMAP).
#'
#' @param n_neighbors Integer. Number of nearest neighbors configuration for the UMAP algorithm.
#' @param min_dist Numeric. Controls how tightly UMAP packs points together in the layout.
#' @param seed Integer. Random seed passed to the UMAP algorithm to ensure reproducible layouts across calls. Default is \code{42}.
#'
#' @export
plotSpectraUMAP <- function(object, n_neighbors = 15, min_dist = 0.1, show_proba = FALSE,
                            palette = NULL, seed = 42, ...) {

  cf               <- umap::umap.defaults
  cf$n_neighbors   <- n_neighbors
  cf$min_dist      <- min_dist
  cf$random_state  <- seed

  set.seed(seed)
  umap_res <- umap::umap(as.matrix(unclass(object$pls_model$scores)),config = cf, ...)
  df <- data.frame(umap_res$layout)
  colnames(df) <- c("U1", "U2")
  df$strain <- as.factor(object$labels)


  lda_alpha <- .computeLDAalpha(object)
  df$Alpha  <- if (show_proba) lda_alpha else rep(0.8, nrow(df))

  n_groups <- nlevels(df$strain)
  pal      <- .getPalette(n_groups, palette)

  range_u1 <- range(df$U1)
  range_u2 <- range(df$U2)
  xlims <- c(range_u1[1] - 0.1 * diff(range_u1), range_u1[2] + 0.1 * diff(range_u1))
  ylims <- c(range_u2[1] - 0.1 * diff(range_u2), range_u2[2] + 0.1 * diff(range_u2))

  p <- ggplot2::ggplot(df, ggplot2::aes(x = U1, y = U2, color = strain)) +
    ggplot2::geom_point(ggplot2::aes(alpha = Alpha), size = 2) +
    ggplot2::scale_alpha_identity() +
    ggplot2::theme_light() +
    ggplot2::scale_color_manual(values = pal) +
    ggplot2::labs(title    = "UMAP projection of spectral profiles",
      subtitle = sprintf("Parameters: n_neighbors = %d, min_dist = %.2f", n_neighbors, min_dist), x = "UMAP 1", y = "UMAP 2") +
    ggplot2::coord_cartesian(xlim = xlims, ylim = ylims)

  return(p)

}
