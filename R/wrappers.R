#' Infer Trajectory with Topological Constraints
#'
#' A high-level wrapper to infer lineages and pseudotime. Calculates centroids
#' from reduced dimension coordinates, constructs a constrained MST, and identifies
#' lineage paths from a root cluster.
#'
#' @param reduced_dim A numeric matrix of coordinates (e.g., PCA, UMAP). Rows are cells, columns are dimensions.
#' @param cluster_labels A vector of cluster IDs for each cell.
#' @param constraints A data.frame of topological constraints (cols: 'from', 'to', 'type').
#' @param start_cluster (Optional) Character string specifying the root cluster.
#' @param probabilistic Logical (default FALSE). If TRUE, uses soft Bayesian priors.
#' @param prior_strength Numeric (default 10). Strength of the prior if probabilistic is TRUE.
#'
#' @return A list containing:
#'   \item{centroids}{Matrix of cluster centroids.}
#'   \item{mst}{The computed Minimum Spanning Tree (igraph object).}
#'   \item{lineages}{A named list of character vectors, where each vector is a sequence of clusters.}
#'   \item{root}{The identifier of the root cluster.}
#'
#' @importFrom igraph degree shortest_paths
#' @export
infer_trajectory <- function(reduced_dim, cluster_labels, constraints = NULL,
                             start_cluster = NULL, probabilistic = FALSE, prior_strength = 10) {

  centroids <- get_centroids(reduced_dim, cluster_labels)
  mst <- get_constrained_mst(centroids, constraints,
                             probabilistic = probabilistic,
                             prior_strength = prior_strength)

  deg <- igraph::degree(mst)
  endpoints <- names(deg)[deg == 1]

  if (!is.null(start_cluster)) {
    root <- start_cluster
    if (!root %in% rownames(centroids)) stop("start_cluster not found in cluster labels.")
  } else {
    root <- if(!is.null(constraints)) as.character(constraints$from[1]) else endpoints[1]
  }

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
