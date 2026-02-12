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
    # robust colMeans for both single and multiple dimensions
    subset_data <- data[labels == cl, , drop = FALSE]
    colMeans(subset_data)
  }))
  
  rownames(centroids) <- unique_clusters
  return(centroids)
}

#' Infer Trajectory with Topological Constraints
#'
#' A high-level wrapper to infer lineages and pseudotime. Calculates centroids
#' from reduced dimension coordinates, constructs a constrained MST using the
#' robust `get_constrained_mst` core, and identifies lineage paths.
#'
#' @param reduced_dim A numeric matrix or data.frame of coordinates (e.g., PCA, UMAP).
#' @param cluster_labels A vector of cluster IDs for each cell.
#' @param constraints A data.frame of topological constraints (cols: 'from', 'to', 'type').
#' @param start_cluster (Optional) Character string specifying the root cluster.
#' @param probabilistic Logical (default FALSE). If TRUE, uses soft Bayesian priors.
#' @param prior_strength Numeric (default 10). Strength of the prior/penalty.
#' @param small_eps Finite positive number for hard 'must_link' (default 1e-8).
#' @param large_value Finite positive number for hard 'cannot_link' (default 1e12).
#'
#' @return A list containing:
#'   \item{centroids}{Matrix of cluster centroids.}
#'   \item{mst}{The computed Minimum Spanning Tree (igraph object).}
#'   \item{lineages}{A named list of character vectors (sequences of clusters).}
#'   \item{root}{The identifier of the root cluster.}
#'
#' @importFrom igraph degree shortest_paths
#' @export
infer_trajectory <- function(reduced_dim, cluster_labels, constraints = NULL, 
                             start_cluster = NULL, 
                             probabilistic = FALSE, prior_strength = 10,
                             small_eps = 1e-8, large_value = 1e12) {
  
  # Calculate Centroids
  centroids <- get_centroids(reduced_dim, cluster_labels)
  
  # Compute Constrained MST (Passing all new parameters)
  mst <- get_constrained_mst(centroids, constraints = constraints, 
                             probabilistic = probabilistic, 
                             prior_strength = prior_strength,
                             small_eps = small_eps,
                             large_value = large_value,
                             return_format = "igraph")
  
  # Identify Endpoints and Root
  deg <- igraph::degree(mst)
  endpoints <- names(deg)[deg == 1]
  
  if (!is.null(start_cluster)) {
    root <- as.character(start_cluster)
    if (!root %in% rownames(centroids)) stop("start_cluster not found in cluster labels.")
  } else {
    # Default: Use first constraint 'from' node or first endpoint
    root <- if(!is.null(constraints)) as.character(constraints$from[1]) else endpoints[1]
  }
  
  # Trace Lineages (Shortest paths on the tree)
  lineages <- list()
  for (end_node in endpoints) {
    if (end_node != root) {
      path_nodes <- igraph::shortest_paths(mst, from = root, to = end_node)$vpath[[1]]
      line_name <- paste0("Lineage_", end_node)
      lineages[[line_name]] <- names(path_nodes)
    }
  }
  
  return(list(
    centroids = centroids,
    mst = mst,
    lineages = lineages,
    root = root
  ))
}
