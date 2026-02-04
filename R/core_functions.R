#' Compute Constrained Minimum Spanning Tree
#'
#' Derives a Minimum Spanning Tree (MST) from cluster centroids while respecting
#' user-defined topological constraints and temporal penalties.
#'
#' @param centroids A numeric matrix of cluster centroids (rows = clusters, cols = dimensions).
#' @param constraints A data.frame with columns 'from', 'to', and 'type'.
#'   Type must be either "must_link" or "cannot_link".
#' @param time_labels A named numeric vector indicating the time point for each cluster.
#'   Used to penalize edges that move backwards in time.
#'
#' @return An \code{igraph} object representing the constrained MST.
#'
#' @importFrom igraph graph_from_adjacency_matrix mst
#' @importFrom stats dist
#' @export
get_constrained_mst <- function(centroids, constraints = NULL, time_labels = NULL) {

  # Calculate base Euclidean distances
  dist_mat <- as.matrix(stats::dist(centroids))
  cluster_names <- rownames(dist_mat)

  # Apply temporal penalties (soft constraints)
  if (!is.null(time_labels)) {
    # Vectorized check for time consistency could replace loops for speed in C++,
    # but R loop is sufficient for cluster-level (N < 100) operations.
    for (i in cluster_names) {
      for (j in cluster_names) {
        if (i == j) next

        # Penalize backward transitions
        if (!is.na(time_labels[i]) && !is.na(time_labels[j])) {
          if (time_labels[j] < time_labels[i]) {
            dist_mat[i, j] <- dist_mat[i, j] * 10
            dist_mat[j, i] <- dist_mat[j, i] * 10
          }
        }
      }
    }
  }

  # Apply topological constraints (hard constraints)
  if (!is.null(constraints)) {
    # Ensure constraint clusters exist in data
    valid_cons <- constraints[constraints$from %in% cluster_names &
                                constraints$to %in% cluster_names, ]

    if (nrow(valid_cons) > 0) {
      must_link <- valid_cons[valid_cons$type == "must_link", ]
      cannot_link <- valid_cons[valid_cons$type == "cannot_link", ]

      # Enforce connections by minimizing distance
      if (nrow(must_link) > 0) {
        dist_mat[cbind(must_link$from, must_link$to)] <- 1e-6
        dist_mat[cbind(must_link$to, must_link$from)] <- 1e-6
      }

      # Sever connections by maximizing distance
      if (nrow(cannot_link) > 0) {
        dist_mat[cbind(cannot_link$from, cannot_link$to)] <- Inf
        dist_mat[cbind(cannot_link$to, cannot_link$from)] <- Inf
      }
    }
  }

  g <- igraph::graph_from_adjacency_matrix(dist_mat, mode = "undirected", weighted = TRUE)
  mst_graph <- igraph::mst(g)

  return(mst_graph)
}
