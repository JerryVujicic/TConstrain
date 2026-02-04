#' Compute Constrained Minimum Spanning Tree
#'
#' Derives a Minimum Spanning Tree (MST) from cluster centroids while respecting
#' user-defined topological constraints and temporal penalties. Supports both
#' hard constraints and soft Bayesian priors.
#'
#' @param centroids A numeric matrix of cluster centroids (rows = clusters, cols = dimensions).
#' @param constraints A data.frame with columns 'from', 'to', and 'type'.
#'   Type must be either "must_link" or "cannot_link".
#' @param time_labels A named numeric vector indicating the time point for each cluster.
#' @param probabilistic Logical. If TRUE, treats constraints as soft priors scaled by
#'   \code{prior_strength}. If FALSE, enforces constraints with effectively zero or infinite distance.
#' @param prior_strength Numeric (default 10). The scaling factor for soft priors.
#'
#' @return An \code{igraph} object representing the constrained MST.
#'
#' @importFrom igraph graph_from_adjacency_matrix mst
#' @importFrom stats dist
#' @export
get_constrained_mst <- function(centroids, constraints = NULL, time_labels = NULL,
                                probabilistic = FALSE, prior_strength = 10) {

  dist_mat <- as.matrix(stats::dist(centroids))
  cluster_names <- rownames(dist_mat)

  if (!is.null(time_labels)) {
    for (i in cluster_names) {
      for (j in cluster_names) {
        if (i == j) next

        if (!is.na(time_labels[i]) && !is.na(time_labels[j]) && time_labels[j] < time_labels[i]) {
          scale_factor <- if (probabilistic) prior_strength else 10
          dist_mat[i, j] <- dist_mat[i, j] * scale_factor
          dist_mat[j, i] <- dist_mat[j, i] * scale_factor
        }
      }
    }
  }

  if (!is.null(constraints)) {
    valid_cons <- constraints[constraints$from %in% cluster_names &
                                constraints$to %in% cluster_names, ]

    if (nrow(valid_cons) > 0) {
      for (k in 1:nrow(valid_cons)) {
        u <- as.character(valid_cons$from[k])
        v <- as.character(valid_cons$to[k])
        ctype <- valid_cons$type[k]

        if (probabilistic) {
          if (ctype == "must_link") {
            new_dist <- dist_mat[u, v] / prior_strength
            dist_mat[u, v] <- new_dist
            dist_mat[v, u] <- new_dist
          } else if (ctype == "cannot_link") {
            new_dist <- dist_mat[u, v] * prior_strength
            dist_mat[u, v] <- new_dist
            dist_mat[v, u] <- new_dist
          }
        } else {
          if (ctype == "must_link") {
            dist_mat[u, v] <- 1e-6
            dist_mat[v, u] <- 1e-6
          } else if (ctype == "cannot_link") {
            dist_mat[u, v] <- Inf
            dist_mat[v, u] <- Inf
          }
        }
      }
    }
  }

  g <- igraph::graph_from_adjacency_matrix(dist_mat, mode = "undirected", weighted = TRUE)
  mst_graph <- igraph::mst(g)

  return(mst_graph)
}
