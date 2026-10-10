#' Get cluster centroids from single cell coordinates
#'
#' @param data Numeric finite matrix or data.frame, cells by dimensions.
#' @param labels One nonmissing, nonempty cluster label per cell.
#' @return A matrix with natural ordering of numeric cluster suffixes. A single
#'   dimension or a single cluster remains a matrix.
#' @export
get_centroids <- function(data, labels) {
  data <- .tc_coordinates(data, named = FALSE)
  if (length(labels) != nrow(data)) stop("'labels' must contain one value for each row of 'data'.")
  labels <- as.character(labels)
  if (anyNA(labels) || any(!nzchar(trimws(labels)))) stop("Labels must be nonmissing and nonempty.")
  clusters <- unique(labels)
  has_suffix <- grepl("[0-9]+$", clusters)
  prefix <- sub("[0-9]+$", "", clusters)
  suffix <- suppressWarnings(as.numeric(sub(".*?([0-9]+)$", "\\1", clusters, perl = TRUE)))
  suffix[!has_suffix] <- 0
  clusters <- clusters[order(prefix, has_suffix, suffix, clusters)]
  centroids <- do.call(rbind, lapply(clusters, function(cl) colMeans(data[labels == cl, , drop = FALSE])))
  rownames(centroids) <- clusters
  centroids
}

#' Infer a rooted cluster topology with exact or probabilistic constraints
#'
#' @param reduced_dim Numeric finite cells by dimensions coordinates.
#' @param cluster_labels One cluster ID per cell.
#' @param start_cluster Explicit root ID. Otherwise use a unique earliest known
#'   time. A single cluster is its own root.
#' @inheritParams get_constrained_mst
#' @return List with centroids, mst, lineages, root, audit and posterior (NULL in
#'   hard mode). Lineages are cluster paths, not cell pseudotimes. Paths are
#'   extracted without optimization weights, which may be negative log scores.
#' @export
infer_trajectory <- function(reduced_dim, cluster_labels, constraints = NULL,
                             time_labels = NULL, start_cluster = NULL,
                             probabilistic = FALSE, prior_strength = 10,
                             small_eps = 1e-8, large_value = 1e12,
                             likelihood_scale = 1, time_strength = 0,
                             time_scale = 1, hard_constraints = NULL) {
  centroids <- get_centroids(reduced_dim, cluster_labels)
  ids <- rownames(centroids)
  times <- .tc_times(time_labels, ids)
  if (!is.null(start_cluster)) {
    root <- as.character(start_cluster)
  } else if (length(ids) == 1L) {
    root <- ids
  } else {
    known <- which(!is.na(times))
    earliest <- if (length(known)) ids[known[times[known] == min(times[known])]] else character()
    if (length(earliest) != 1L) stop("Provide 'start_cluster' or time_labels with a unique earliest known time.")
    root <- earliest
  }
  if (length(root) != 1L || is.na(root) || !root %in% ids) stop("'start_cluster' must name one inferred cluster.")
  args <- list(centroids = centroids, constraints = constraints, time_labels = time_labels,
               probabilistic = probabilistic, prior_strength = prior_strength,
               likelihood_scale = likelihood_scale, time_strength = time_strength,
               time_scale = time_scale, hard_constraints = hard_constraints)
  if (!missing(small_eps)) args$small_eps <- small_eps
  if (!missing(large_value)) args$large_value <- large_value
  mst <- do.call(get_constrained_mst, args)
  if (length(ids) == 1L) {
    lineages <- list(Lineage_1 = root)
  } else {
    degree <- igraph::degree(mst)
    endpoints <- setdiff(names(degree)[degree == 1L], root)
    lineages <- lapply(endpoints, function(end) {
      names(igraph::shortest_paths(mst, from = root, to = end, weights = NA)$vpath[[1]])
    })
    names(lineages) <- paste0("Lineage_", endpoints)
  }
  list(centroids = centroids, mst = mst, lineages = lineages, root = root,
       audit = igraph::graph_attr(mst, "audit"), posterior = igraph::graph_attr(mst, "posterior"))
}
