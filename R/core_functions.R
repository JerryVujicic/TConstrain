#' Compute Constrained Minimum Spanning Tree (validated, configurable, tidy-friendly)
#'
#' Builds a Minimum Spanning Tree (MST) on cluster centroids while honoring
#' temporal directionality and user-provided topological constraints. Constraints
#' can be applied as hard rules or as soft (probabilistic) priors.
#'
#' @param centroids A numeric matrix or data.frame/tibble of cluster centroids
#'   (rows = clusters, columns = dimensions). **Row names must be cluster IDs**
#'   unless using the tidy wrapper which can accept an explicit `id_col`.
#' @param constraints Optional data.frame (or tibble) with columns: `from`,
#'   `to`, and `type`. `type` must be one of `"must_link"` or `"cannot_link"`.
#' @param time_labels Optional **named** numeric vector. Names must match the
#'   row names of `centroids`. If provided, edges that go "back in time"
#'   (from a cluster at time t1 to a cluster at an earlier time t0) will be
#'   penalized. NA times are treated as "unknown" and do not trigger penalties.
#' @param probabilistic Logical scalar. If TRUE, constraints and temporal
#'   penalties are applied as *soft priors* (distances are scaled by
#'   `prior_strength`). If FALSE, constraints are enforced as hard rules
#'   using `small_eps` and `large_value`.
#' @param prior_strength Numeric > 1. Strength of prior/penalty used when
#'   `probabilistic = TRUE`; also used as the base scale for the non-probabilistic
#'   temporal penalty (for consistent behaviour).
#' @param small_eps Finite positive number used as the effective distance for a
#'   hard `must_link`. Defaults to `1e-8`.
#' @param large_value Finite positive number used as the effective distance for a
#'   hard `cannot_link`. Defaults to `1e12`.
#' @param return_format Character: either `"igraph"` (default) or `"edgelist"`.
#' @return Either an `igraph` object or an edge list data.frame with columns
#'   `from`, `to`, and `weight`.
#' @seealso get_constrained_mst_tidy
#' @importFrom igraph graph_from_adjacency_matrix mst E V
#' @importFrom stats dist
#' @export
get_constrained_mst <- function(centroids,
                                constraints = NULL,
                                time_labels = NULL,
                                probabilistic = FALSE,
                                prior_strength = 10,
                                small_eps = 1e-8,
                                large_value = 1e12,
                                return_format = c("igraph", "edgelist")) {
  return_format <- match.arg(return_format)

  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' is required but not installed. Please install it.")
  }
  if (!is.matrix(centroids) && !is.data.frame(centroids)) {
    stop("'centroids' must be a numeric matrix or data.frame with row names (cluster IDs).")
  }
  centroids <- as.matrix(centroids)
  if (!is.numeric(centroids)) stop("'centroids' must be numeric.")
  cluster_names <- rownames(centroids)
  if (is.null(cluster_names) || any(cluster_names == "")) {
    stop("Row names of 'centroids' must be non-empty cluster identifiers.")
  }
  if (!is.logical(probabilistic) || length(probabilistic) != 1) stop("'probabilistic' must be a single logical value (TRUE/FALSE).")
  if (!is.numeric(prior_strength) || length(prior_strength) != 1 || prior_strength <= 1) stop("'prior_strength' must be a single numeric value > 1.")
  if (!is.numeric(small_eps) || length(small_eps) != 1 || small_eps <= 0) stop("'small_eps' must be a positive numeric scalar.")
  if (!is.numeric(large_value) || length(large_value) != 1 || large_value <= 0) stop("'large_value' must be a positive numeric scalar.")

  dist_mat <- as.matrix(stats::dist(centroids))
  if (!is.matrix(dist_mat) || nrow(dist_mat) != ncol(dist_mat)) stop("Failed to compute a square distance matrix.")
  n <- nrow(dist_mat)
  colnames(dist_mat) <- rownames(dist_mat) <- cluster_names

  if (!is.null(time_labels)) {
    if (!is.numeric(time_labels) || is.null(names(time_labels))) {
      stop("'time_labels' must be a named numeric vector with names matching cluster IDs (rownames of centroids).")
    }
    unknown_time_names <- setdiff(names(time_labels), cluster_names)
    if (length(unknown_time_names) > 0) warning(sprintf("Ignoring %d time labels with unknown cluster IDs.", length(unknown_time_names)))
    times_vec <- time_labels[cluster_names]
    time_back <- outer(times_vec, times_vec, FUN = function(ti, tj) ifelse(is.na(ti) | is.na(tj), FALSE, tj < ti))
    scale_factor <- if (probabilistic) prior_strength else max(prior_strength, 10)
    multiplier <- matrix(1, nrow = n, ncol = n)
    multiplier[time_back] <- scale_factor
    diag(multiplier) <- 1
    dist_mat <- dist_mat * multiplier
  }

  if (!is.null(constraints)) {
    if (!is.data.frame(constraints)) stop("'constraints' must be a data.frame with columns 'from', 'to', 'type'.")
    req_cols <- c("from", "to", "type")
    if (!all(req_cols %in% colnames(constraints))) stop(paste0("'constraints' must contain columns: ", paste(req_cols, collapse = ", ")))
    constraints$from <- as.character(constraints$from)
    constraints$to <- as.character(constraints$to)
    constraints$type <- as.character(constraints$type)
    keep_idx <- which(constraints$from %in% cluster_names & constraints$to %in% cluster_names)
    if (length(keep_idx) < nrow(constraints)) warning(sprintf("Dropping %d constraints with unknown cluster IDs.", nrow(constraints) - length(keep_idx)))
    valid_cons <- constraints[keep_idx, , drop = FALSE]
    allowed_types <- c("must_link", "cannot_link")
    bad_types <- setdiff(unique(valid_cons$type), allowed_types)
    if (length(bad_types) > 0) stop(sprintf("Unknown constraint types found: %s. Allowed: %s", paste(bad_types, collapse = ", "), paste(allowed_types, collapse = ", ")))
    if (nrow(valid_cons) > 0) {
      pair_key <- paste(pmin(valid_cons$from, valid_cons$to), pmax(valid_cons$from, valid_cons$to), sep = "--")
      duplicate_pairs <- unique(pair_key[duplicated(pair_key)])
      if (length(duplicate_pairs) > 0) warning(sprintf("Duplicate constraint pairs found; later rows override earlier rows: %s", paste(duplicate_pairs, collapse = ", ")))
      for (k in seq_len(nrow(valid_cons))) {
        u <- valid_cons$from[k]
        v <- valid_cons$to[k]
        ctype <- valid_cons$type[k]
        if (probabilistic) {
          if (ctype == "must_link") {
            dist_mat[u, v] <- dist_mat[u, v] / prior_strength
            dist_mat[v, u] <- dist_mat[v, u] / prior_strength
          } else if (ctype == "cannot_link") {
            dist_mat[u, v] <- dist_mat[u, v] * prior_strength
            dist_mat[v, u] <- dist_mat[v, u] * prior_strength
          }
        } else if (ctype == "must_link") {
          dist_mat[u, v] <- small_eps
          dist_mat[v, u] <- small_eps
        } else if (ctype == "cannot_link") {
          dist_mat[u, v] <- large_value
          dist_mat[v, u] <- large_value
        }
      }
    }
  }

  diag(dist_mat) <- 0
  is_bad <- !is.finite(dist_mat)
  if (any(is_bad)) {
    dist_mat[is_bad] <- large_value
    warning("Some distances were non-finite (NA/Inf) and were replaced by 'large_value'. This may disconnect the graph.")
  }
  dist_mat <- (dist_mat + t(dist_mat)) / 2
  g <- igraph::graph_from_adjacency_matrix(dist_mat, mode = "undirected", weighted = TRUE, diag = FALSE)
  mst_graph <- igraph::mst(g)
  igraph::graph_attr(mst_graph, "modified_distance_matrix") <- dist_mat
  if (return_format == "igraph") return(mst_graph)
  edges <- igraph::as_data_frame(mst_graph, what = "edges")
  edges <- edges[, c("from", "to", "weight"), drop = FALSE]
  return(edges)
}

#' Tidy wrapper: accept a tibble/data.frame of centroids (with an ID column)
#'
#' @param data A data.frame/tibble of centroid coordinates + an ID column.
#' @param id_col Character. Column name containing cluster IDs (will become rownames).
#' @inheritParams get_constrained_mst
#' @export
get_constrained_mst_tidy <- function(data, id_col = "id", constraints = NULL,
                                     time_labels = NULL, probabilistic = FALSE,
                                     prior_strength = 10, small_eps = 1e-8,
                                     large_value = 1e12, return_format = c("igraph", "edgelist")) {
  return_format <- match.arg(return_format)
  if (!is.data.frame(data)) stop("'data' must be a data.frame or tibble.")
  if (!id_col %in% colnames(data)) stop(sprintf("id_col '%s' not found in data.", id_col))
  cent <- as.data.frame(data)
  rownames(cent) <- as.character(cent[[id_col]])
  cent[[id_col]] <- NULL
  get_constrained_mst(centroids = cent, constraints = constraints,
                      time_labels = time_labels, probabilistic = probabilistic,
                      prior_strength = prior_strength, small_eps = small_eps,
                      large_value = large_value, return_format = return_format)
}

#' Run simple simulated examples and checks to validate behaviour
#'
#' @param verbose Logical; if TRUE prints intermediate results.
#' @return TRUE invisibly when checks pass; otherwise stops with an error.
#' @export
run_constrained_mst_checks <- function(verbose = TRUE) {
  if (verbose) message("Running simple simulated checks...")
  set.seed(42)
  cent <- matrix(rnorm(10), nrow = 5)
  rownames(cent) <- paste0("C", 1:5)
  times <- c(C1 = 0, C2 = 1, C3 = 2, C4 = 1, C5 = 3)
  cons <- data.frame(from = c("C1", "C4"), to = c("C2", "C5"), type = c("must_link", "cannot_link"), stringsAsFactors = FALSE)
  mst_hard <- get_constrained_mst(cent, constraints = cons, time_labels = times, probabilistic = FALSE, prior_strength = 5, return_format = "edgelist")
  if (verbose) print(mst_hard)
  e_names <- apply(mst_hard[, c("from", "to")], 1, function(x) paste(sort(x), collapse = "--"))
  if (!any(e_names == "C1--C2")) stop("Hard must_link failed: edge C1--C2 not present in MST.")
  mst_soft_graph <- get_constrained_mst(cent, constraints = cons, time_labels = times, probabilistic = TRUE, prior_strength = 10, return_format = "igraph")
  if (!inherits(mst_soft_graph, "igraph")) stop("Expected igraph when return_format = 'igraph'.")
  if (verbose) message("All checks passed.")
  invisible(TRUE)
}
