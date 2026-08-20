#' Get Cluster Centroids from Single-Cell Data
#'
#' Calculates the mean coordinate of each cluster. Labels with numeric suffixes
#' use natural ordering, so `C2` precedes `C10`.
#'
#' @param data A numeric matrix or data.frame of coordinates (cells x dims).
#' @param labels A vector of cluster labels corresponding to cells.
#' @return A matrix of cluster centroids (clusters x dims).
#' @export
get_centroids <- function(data, labels) {
  data <- as.matrix(data)
  if (nrow(data) != length(labels)) {
    stop("'labels' must contain one value for each row of 'data'.")
  }
  labels <- as.character(labels)
  unique_clusters <- unique(labels)
  has_numeric_suffix <- grepl("[0-9]+$", unique_clusters)
  prefix <- sub("[0-9]+$", "", unique_clusters)
  numeric_suffix <- suppressWarnings(as.numeric(
    sub(".*?([0-9]+)$", "\\1", unique_clusters, perl = TRUE)
  ))
  numeric_suffix[!has_numeric_suffix] <- 0
  unique_clusters <- unique_clusters[order(
    prefix,
    has_numeric_suffix,
    numeric_suffix,
    unique_clusters,
    na.last = TRUE
  )]
  centroids <- t(sapply(unique_clusters, function(cl) {
    subset_data <- data[labels == cl, , drop = FALSE]
    colMeans(subset_data)
  }))
  rownames(centroids) <- unique_clusters
  return(centroids)
}

#' Infer Trajectory with Topological Constraints
#'
#' A high-level wrapper that calculates cluster centroids, constructs a
#' constrained MST with optional temporal penalties, and identifies lineage paths
#' from a root cluster. Constraint-table order is never used to choose the root.
#'
#' @param reduced_dim A numeric matrix or data.frame of coordinates (e.g., PCA, UMAP).
#' @param cluster_labels A vector of cluster IDs for each cell.
#' @param constraints A data.frame of topological constraints (cols: 'from', 'to', 'type').
#' @param time_labels Optional named numeric vector assigning a time point to each cluster.
#' @param start_cluster Optional character string specifying the root cluster.
#'   When omitted, a unique cluster with the earliest known `time_labels` value
#'   is used. It must be supplied when no unique earliest time label is available.
#' @param probabilistic Logical (default FALSE). If TRUE, uses soft Bayesian priors.
#' @param prior_strength Numeric (default 10). Strength of the prior/penalty.
#' @param small_eps Finite positive number for hard 'must_link' (default 1e-8).
#' @param large_value Finite positive number for hard 'cannot_link' (default 1e12).
#' @return A list containing centroids, MST, lineages, and root.
#' @export
infer_trajectory <- function(reduced_dim, cluster_labels, constraints = NULL,
                             time_labels = NULL, start_cluster = NULL,
                             probabilistic = FALSE, prior_strength = 10,
                             small_eps = 1e-8, large_value = 1e12) {
  centroids <- get_centroids(reduced_dim, cluster_labels)

  if (!is.null(start_cluster)) {
    root <- as.character(start_cluster)
  } else {
    if (is.null(time_labels)) {
      stop("Provide 'start_cluster' or named 'time_labels' with a unique earliest time to infer the trajectory root.")
    }
    if (!is.numeric(time_labels) || is.null(names(time_labels))) {
      stop("'time_labels' must be a named numeric vector when 'start_cluster' is not supplied.")
    }
    aligned_times <- time_labels[rownames(centroids)]
    known_times <- aligned_times[!is.na(aligned_times)]
    if (!length(known_times)) {
      stop("Cannot infer a root because 'time_labels' has no known values for the supplied clusters. Provide 'start_cluster'.")
    }
    earliest_clusters <- names(known_times)[known_times == min(known_times)]
    if (length(earliest_clusters) != 1) {
      stop("Cannot infer a unique root from 'time_labels'. Provide 'start_cluster'.")
    }
    root <- earliest_clusters
  }

  if (length(root) != 1 || is.na(root) || !root %in% rownames(centroids)) {
    stop("'start_cluster' must name one of the inferred clusters.")
  }

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

  lineages <- list()
  for (end_node in endpoints) {
    if (end_node != root) {
      path_nodes <- igraph::shortest_paths(mst, from = root, to = end_node)$vpath[[1]]
      lineages[[paste0("Lineage_", end_node)]] <- names(path_nodes)
    }
  }

  list(centroids = centroids, mst = mst, lineages = lineages, root = root)
}
