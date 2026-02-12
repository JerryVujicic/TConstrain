#' Compute Constrained Minimum Spanning Tree (validated, configurable, tidy-friendly)
#'
#' Builds a Minimum Spanning Tree (MST) on cluster centroids while honoring
#' temporal directionality and user-provided topological constraints. Constraints
#' can be applied as hard rules or as soft (probabilistic) priors.
#'
#' This version adds:
#' * explicit parameters for hard/soft magnitudes (`small_eps`, `large_value`),
#' * a `return_format` option to get either an `igraph` object or a tidy edge list,
#' * a tidyverse-friendly wrapper that accepts tibbles and returns tibbles,
#' * built-in simulated examples and simple unit checks (no external testthat required),
#' * more comprehensive documentation and friendly messages.
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
#' @param return_format Character: either `"igraph"` (default) to return an
#'   `igraph` object, or `"edgelist"` to return a tidy tibble with columns
#'   `from`, `to`, and `weight` describing the MST edges.
#'
#' @return Either an `igraph` object (when `return_format = "igraph"`) or a
#'   tibble edge list (`return_format = "edgelist"`). When returning an igraph,
#'   the modified distance matrix is stored as a graph attribute
#'   `modified_distance_matrix` for inspection.
#'
#' @seealso get_constrained_mst_tidy for a wrapper that accepts tibbles and an
#'   explicit `id_col`.
#'
#' @examples
#' # Example (run manually):
#' # Basic usage
#' cent <- matrix(rnorm(30), nrow = 6)
#' rownames(cent) <- paste0("C", 1:6)
#' times <- c(C1 = 0, C2 = 1, C3 = 2, C4 = 1, C5 = 3, C6 = 2)
#' cons <- data.frame(from = c("C1","C3"), to = c("C2","C6"), type = c("must_link","cannot_link"))
#' mst_ig <- get_constrained_mst(cent, constraints = cons, time_labels = times,
#'                               probabilistic = TRUE, prior_strength = 5)
#' # Get tidy edge list instead
#' mst_tbl <- get_constrained_mst(cent, constraints = cons, time_labels = times,
#'                                probabilistic = FALSE, return_format = "edgelist")
#'
#' # Tidy wrapper (data.frame/tibble input)
#' library(tibble)
#' cent_tbl <- as_tibble(cent)
#' cent_tbl$id <- rownames(cent)
#' mst_from_tidy <- get_constrained_mst_tidy(cent_tbl, id_col = "id", constraints = cons,
#'                                           time_labels = times, return_format = "edgelist")
#'
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

  ## Basic validation
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

  if (!is.logical(probabilistic) || length(probabilistic) != 1) {
    stop("'probabilistic' must be a single logical value (TRUE/FALSE).")
  }
  if (!is.numeric(prior_strength) || length(prior_strength) != 1 || prior_strength <= 1) {
    stop("'prior_strength' must be a single numeric value > 1.")
  }
  if (!is.numeric(small_eps) || length(small_eps) != 1 || small_eps <= 0) stop("'small_eps' must be a positive numeric scalar.")
  if (!is.numeric(large_value) || length(large_value) != 1 || large_value <= 0) stop("'large_value' must be a positive numeric scalar.")

  ## Compute pairwise distances (Euclidean)
  dist_mat <- as.matrix(stats::dist(centroids))
  if (!is.matrix(dist_mat) || nrow(dist_mat) != ncol(dist_mat)) stop("Failed to compute a square distance matrix.")
  n <- nrow(dist_mat)
  colnames(dist_mat) <- rownames(dist_mat) <- cluster_names

  ## Temporal penalties (vectorized)
  if (!is.null(time_labels)) {
    if (!is.numeric(time_labels)) stop("'time_labels' must be a named numeric vector.")
    if (is.null(names(time_labels))) stop("'time_labels' must be a named numeric vector with names matching cluster IDs (rownames of centroids).")

    times_vec <- time_labels[cluster_names]

    time_back <- outer(times_vec, times_vec, FUN = function(ti, tj) {
      # if either is NA produce FALSE (no penalty)
      ifelse(is.na(ti) | is.na(tj), FALSE, tj < ti)
    })

    scale_factor <- if (probabilistic) prior_strength else max(prior_strength, 10)

    multiplier <- matrix(1, nrow = n, ncol = n)
    multiplier[time_back] <- scale_factor
    diag(multiplier) <- 1

    dist_mat <- dist_mat * multiplier
  }

  ## Constraints handling
  if (!is.null(constraints)) {
    if (!is.data.frame(constraints)) stop("'constraints' must be a data.frame with columns 'from', 'to', 'type'.")
    req_cols <- c("from", "to", "type")
    if (!all(req_cols %in% colnames(constraints))) stop(paste0("'constraints' must contain columns: ", paste(req_cols, collapse = ", ")))

    constraints$from <- as.character(constraints$from)
    constraints$to   <- as.character(constraints$to)
    constraints$type <- as.character(constraints$type)

    keep_idx <- which(constraints$from %in% cluster_names & constraints$to %in% cluster_names)
    if (length(keep_idx) < nrow(constraints)) {
      warning(sprintf("Dropping %d constraints with unknown cluster IDs.", nrow(constraints) - length(keep_idx)))
    }
    valid_cons <- constraints[keep_idx, , drop = FALSE]

    allowed_types <- c("must_link", "cannot_link")
    bad_types <- setdiff(unique(valid_cons$type), allowed_types)
    if (length(bad_types) > 0) stop(sprintf("Unknown constraint types found: %s. Allowed: %s",
                                           paste(bad_types, collapse = ", "), paste(allowed_types, collapse = ", ")))

    if (nrow(valid_cons) > 0) {
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
        } else {
          if (ctype == "must_link") {
            dist_mat[u, v] <- small_eps
            dist_mat[v, u] <- small_eps
          } else if (ctype == "cannot_link") {
            dist_mat[u, v] <- large_value
            dist_mat[v, u] <- large_value
          }
        }
      }
    }
  }

  ## Final sanitization before graph construction
  diag(dist_mat) <- 0
  is_bad <- !is.finite(dist_mat)
  if (any(is_bad)) {
    dist_mat[is_bad] <- large_value
    warning("Some distances were non-finite (NA/Inf) and were replaced by 'large_value'. This may disconnect the graph.")
  }
  dist_mat <- (dist_mat + t(dist_mat)) / 2

  ## Build igraph and compute MST
  g <- igraph::graph_from_adjacency_matrix(dist_mat, mode = "undirected", weighted = TRUE, diag = FALSE)
  mst_graph <- igraph::mst(g)

  igraph::graph_attr(mst_graph, "modified_distance_matrix") <- dist_mat

  if (return_format == "igraph") return(mst_graph)

  ## Convert to tidy edgelist (if requested)
  edges <- igraph::as_data_frame(mst_graph, what = "edges")
  # as_data_frame returns columns: from (vertex id), to (vertex id), weight
  # map numeric vertex ids to names
  vnames <- igraph::V(mst_graph)$name
  edges$from <- vnames[as.integer(edges$from)]
  edges$to   <- vnames[as.integer(edges$to)]
  # ensure tidy column order
  edges <- edges[, c("from", "to", "weight")]
  # return as data.frame/tibble
  return(edges)
}


#' Tidy wrapper: accept a tibble/data.frame of centroids (with an ID column)
#'
#' This wrapper makes it convenient to use data.frame/tibble centroids where
#' the cluster identifier is stored in a column (e.g., `id_col = "cluster_id").
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

  # Call main function
  out <- get_constrained_mst(centroids = cent, constraints = constraints,
                             time_labels = time_labels, probabilistic = probabilistic,
                             prior_strength = prior_strength, small_eps = small_eps,
                             large_value = large_value, return_format = return_format)
  return(out)
}


#' Run simple simulated examples and checks to validate behaviour
#'
#' This function runs a short self-check: constructs a small simulated centroid
#' set, applies a couple of constraints and time labels, and verifies basic
#' expectations (e.g., must_link produces a tiny-distance connection when
#' non-probabilistic, return formats are correct).
#'
#' @param verbose Logical; if TRUE prints intermediate results.
#' @return TRUE invisibly when checks pass; otherwise stops with an error.
#' @export
run_constrained_mst_checks <- function(verbose = TRUE) {
  if (verbose) message("Running simple simulated checks...")

  # simulate 5 centroids in 2D
  set.seed(42)
  cent <- matrix(rnorm(10), nrow = 5)
  rownames(cent) <- paste0("C", 1:5)

  # simple time labels that encourage edges C3->C1 to be penalized
  times <- c(C1 = 0, C2 = 1, C3 = 2, C4 = 1, C5 = 3)

  # force a must-link between C1 and C2 and forbid C4-C5
  cons <- data.frame(from = c("C1", "C4"), to = c("C2", "C5"), type = c("must_link", "cannot_link"), stringsAsFactors = FALSE)

  # Hard constraints
  mst_hard <- get_constrained_mst(cent, constraints = cons, time_labels = times, probabilistic = FALSE, prior_strength = 5, return_format = "edgelist")
  if (verbose) print(mst_hard)
  # check that C1--C2 exists as an edge (must_link)
  e_names <- apply(mst_hard[, c("from","to")], 1, function(x) paste(sort(x), collapse = "--"))
  if (!any(e_names == "C1--C2")) stop("Hard must_link failed: edge C1--C2 not present in MST.")

  # Soft constraints: prior_strength should reduce distance but not force connectivity
  mst_soft_graph <- get_constrained_mst(cent, constraints = cons, time_labels = times, probabilistic = TRUE, prior_strength = 10, return_format = "igraph")
  if (!inherits(mst_soft_graph, "igraph")) stop("Expected igraph when return_format = 'igraph'.")

  if (verbose) message("All checks passed.")
  invisible(TRUE)
}
