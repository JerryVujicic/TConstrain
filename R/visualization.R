#' Interactive Visualization of Topological Constraints
#'
#' Visualizes the reduced dimension plot with overlaid user constraints.
#' Useful for verifying "Must Link" and "Cannot Link" logic before inference.
#'
#' @param reduced_dim A matrix of coordinates (rows=cells, cols=dims).
#' @param cluster_labels A vector of cluster labels.
#' @param constraints A data.frame of constraints.
#'
#' @return A plotly object.
#' @importFrom ggplot2 ggplot geom_point aes theme_minimal labs geom_segment arrow unit
#' @importFrom plotly ggplotly
#' @export
plot_constraints_interactive <- function(reduced_dim, cluster_labels, constraints = NULL) {

  df <- data.frame(
    Dim1 = reduced_dim[,1],
    Dim2 = reduced_dim[,2],
    Cluster = as.factor(cluster_labels)
  )

  centroids <- get_centroids(reduced_dim, cluster_labels)
  cent_df <- as.data.frame(centroids)
  cent_df$Cluster <- rownames(centroids)

  p <- ggplot2::ggplot() +
    ggplot2::geom_point(data = df, ggplot2::aes(x = Dim1, y = Dim2, color = Cluster),
                        alpha = 0.5, size = 1) +
    ggplot2::geom_point(data = cent_df, ggplot2::aes(x = V1, y = V2),
                        size = 5, shape = 21, fill = "white", stroke = 2) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "Trajectory Constraints Map")

  if (!is.null(constraints)) {
    for(i in 1:nrow(constraints)) {
      u <- as.character(constraints$from[i])
      v <- as.character(constraints$to[i])
      type <- constraints$type[i]

      if(u %in% rownames(centroids) && v %in% rownames(centroids)) {

        seg_data <- data.frame(
          x = centroids[u, 1], y = centroids[u, 2],
          xend = centroids[v, 1], yend = centroids[v, 2]
        )

        if (type == "must_link") {
          p <- p + ggplot2::geom_segment(data = seg_data,
                                         ggplot2::aes(x=x, y=y, xend=xend, yend=yend),
                                         arrow = ggplot2::arrow(length = ggplot2::unit(0.03, "npc")),
                                         color = "darkgreen", size = 1.2)
        } else if (type == "cannot_link") {
          p <- p + ggplot2::geom_segment(data = seg_data,
                                         ggplot2::aes(x=x, y=y, xend=xend, yend=yend),
                                         linetype = "dashed", color = "firebrick", size = 1.2, alpha=0.8)
        }
      }
    }
  }

  return(plotly::ggplotly(p))
}
