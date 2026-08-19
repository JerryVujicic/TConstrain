#' Get Cluster Centroids from Single-Cell Data
#'
#' Helper function to calculate the center of mass for each cluster.
#'
#' @param data A numeric matrix or data.frame of coordinates (cells x dims).
#' @param labels A vector of cluster labels corresponding to cells.
#' @return A matrix of cluster centroids (clusters x dims).
#' @export
get_centroids <- function(data, labels) {
  data <- as.matrix(data)
  unique_clusters <- sort(unique(as.character(labels)))
  centroids <- t(sapply(unique_clusters, function(cl) {
    subset_data <- data[labels == cl, , drop = FALSE]
    colMeans(subset_data)
  }))
  rownames(centroids) <- unique_clusters
  return(centroids)
}

#' Infer Trajectory with Topological Constraints
#'
#' @param reduced_dim A numeric matrix or data.frame of coordinates (e.g., PCA, UMAP).
#' @param cluster_labels A vector of cluster IDs for each cell.
#' @param constraints A data.frame of topological constraints (cols: 'from', 'to', 'type').
#' @param time_labels Optional named numeric vector assigning a time point to each cluster.
#' @param start_cluster Optional character string specifying the root cluster.
#' @param probabilistic Logical (default FALSE). If TRUE, uses soft Bayesian priors.
#' @param prior_strength Numeric (default 10). Strength of the prior/penalty.
#' @param small_eps Finite positive number for hard 'must_link' (default 1e-8).
#' @param large_value Finite positive number for hard 'cannot_link' (default 1e12).
#' @return A list containing centroids, MST, lineages, and root.
#' @importFrom igraph degree shortest_paths
#' @export
infer_trajectory <- function(reduced_dim, cluster_labels, constraints = NULL,
                             time_labels = NULL, start_cluster = NULL,
                             probabilistic = FALSE, prior_strength = 10,
                             small_eps = 1e-8, large_value = 1e12) {
  centroids <- get_centroids(reduced_dim, cluster_labels)

  mst <- get_constrained_mst(centroids, constraints = constraints,
                             time_labels = time_labels,
                             probabilistic = probabilistic,
                             prior_strength = prior_strength,
                             small_eps = small_eps,
                             large_value = large_value,
                             return_format = "igraph")

  deg <- igraph::degree(mst)
  endpoints <- names(deg)[deg == 1]
  if (!length(endpoints)) stop("Could not identify an endpoint in the inferred MST.")

  if (!is.null(start_cluster)) {
    root <- as.character(start_cluster)
  } else if (!is.null(constraints) && nrow(constraints) > 0) {
    root <- as.character(constraints$from[1])
  } else {
    root <- endpoints[1]
  }

  if (length(root) != 1 || is.na(root) || !root %in% rownames(centroids)) {
    stop("Inferred root cluster not found in cluster labels. Check constraints$from[1] or provide a valid start_cluster.")
  }

  lineages <- list()
  for (end_node in endpoints) {
    if (end_node != root) {
      path_nodes <- igraph::shortest_paths(mst, from = root, to = end_node)$vpath[[1]]
      lineages[[paste0("Lineage_", end_node)]] <- names(path_nodes)
    }
  }

  list(centroids = centroids, mst = mst, lineages = lineages, root = root)
}
