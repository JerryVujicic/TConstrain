#' Interactive Visualization of Topological Constraints
#'
#' Visualizes the reduced dimension plot with overlaid user constraints.
#' Useful for verifying "Must Link" and "Cannot Link" logic before inference.
#' Accepts both matrix and data.frame/tibble inputs.
#'
#' @param reduced_dim A matrix or data.frame of coordinates (rows=cells, cols=dims).
#' @param cluster_labels A vector of cluster labels.
#' @param constraints A data.frame of constraints (cols: from, to, type).
#' @return A plotly object.
#' @importFrom ggplot2 ggplot geom_point aes theme_minimal labs geom_segment arrow unit
#' @importFrom plotly ggplotly
#' @export
plot_constraints_interactive <- function(reduced_dim, cluster_labels, constraints = NULL) {
  reduced_mat <- as.matrix(reduced_dim)
  if (ncol(reduced_mat) < 2) stop("'reduced_dim' must contain at least two dimensions.")

  df <- data.frame(
    Dim1 = reduced_mat[, 1],
    Dim2 = reduced_mat[, 2],
    Cluster = as.factor(cluster_labels)
  )

  centroids <- get_centroids(reduced_mat, cluster_labels)
  cent_df <- data.frame(
    Dim1 = centroids[, 1],
    Dim2 = centroids[, 2],
    Cluster = rownames(centroids),
    stringsAsFactors = FALSE
  )

  p <- ggplot2::ggplot() +
    ggplot2::geom_point(data = df, ggplot2::aes(x = Dim1, y = Dim2, color = Cluster),
                        alpha = 0.5, size = 1) +
    ggplot2::geom_point(data = cent_df, ggplot2::aes(x = Dim1, y = Dim2),
                        size = 5, shape = 21, fill = "white", stroke = 2) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "Trajectory Constraints Map")

  if (!is.null(constraints)) {
    constraints$from <- as.character(constraints$from)
    constraints$to <- as.character(constraints$to)

    for (i in seq_len(nrow(constraints))) {
      u <- constraints$from[i]
      v <- constraints$to[i]
      type <- constraints$type[i]
      if (u %in% rownames(centroids) && v %in% rownames(centroids)) {
        seg_data <- data.frame(
          x = centroids[u, 1], y = centroids[u, 2],
          xend = centroids[v, 1], yend = centroids[v, 2]
        )
        if (type == "must_link") {
          p <- p + ggplot2::geom_segment(
            data = seg_data,
            ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
            arrow = ggplot2::arrow(length = ggplot2::unit(0.03, "npc")),
            color = "darkgreen", linewidth = 1.2
          )
        } else if (type == "cannot_link") {
          p <- p + ggplot2::geom_segment(
            data = seg_data,
            ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
            linetype = "dashed", color = "firebrick", linewidth = 1.2, alpha = 0.8
          )
        }
      }
    }
  }

  plotly::ggplotly(p)
}
